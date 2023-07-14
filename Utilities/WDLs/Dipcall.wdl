version 1.0

# Modified from "official" Dockstore version to allow custom PAR regions for male samples
# Allows for use with non-hg37/38 references
workflow runDipcall {
    input {
        File assemblyFastaPat
        File assemblyFastaMat
        File referenceFasta
        Boolean isMaleSample
        File? custom_PAR_bed
        String sample_name = "dipcall"
    }

	call dipcall {
        input:
            assemblyFastaPat=assemblyFastaPat,
            assemblyFastaMat=assemblyFastaMat,
            referenceFasta=referenceFasta,
            isMaleSample=isMaleSample,
            custom_PAR_bed=custom_PAR_bed
	}

    # Note: Highly recommended to include custom PAR bed to clean missing genotypes in haploid regions for male samples
    # This allows easier downstream analysis for comparison tools
    if ((isMaleSample) && (defined(custom_PAR_bed))) {
        call CleanHaploidCalls {
            input:
                dipcall_vcf=dipcall.outputVCF,
                PAR_bed=select_first([custom_PAR_bed])
        }
    }

    call IndexVCF {
        input:
            input_vcf=select_first([CleanHaploidCalls.cleaned_vcf, dipcall.outputVCF]),
            sample_name=sample_name
    }

    output{
        File outputTarball    = dipcall.outputTarball
        File outputVCF        = IndexVCF.vcf
        File outputVCF_index = IndexVCF.vcf_index
        File outputBED        = dipcall.outputBED
    }
}

task dipcall {
    input {
        File assemblyFastaPat
        File assemblyFastaMat
        File referenceFasta
        Boolean isMaleSample
        Boolean referenceIsHS38 = true
        File? custom_PAR_bed
        Int memSizeGB = 64
        Int threadCount = 16
        Int diskSizeGB = 64
        String dockerImage = "humanpangenomics/hpp_dipcall_v0.3:latest"
    }

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace
        PATH="/root/bin/samtools_1.9:$PATH"

        # get output base
        PREFIX=$(basename ~{assemblyFastaPat} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//')
        mkdir $PREFIX.dipcall

        # prep paternal
        PAT_FILENAME=$(basename -- "~{assemblyFastaPat}")
        if [[ $PAT_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFastaPat} .
            gunzip $PAT_FILENAME
            PAT_FILENAME="${PAT_FILENAME%.gz}"
        else
            ln -s ~{assemblyFastaPat}
        fi

        # prep maternal
        MAT_FILENAME=$(basename -- "~{assemblyFastaMat}")
        if [[ $MAT_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFastaMat} .
            gunzip $MAT_FILENAME
            MAT_FILENAME="${MAT_FILENAME%.gz}"
        else
            ln -s ~{assemblyFastaMat}
        fi

        # prep reference
        REF_FILENAME=$(basename -- "~{referenceFasta}")
        if [[ $REF_FILENAME =~ \.gz$ ]]; then
            cp ~{referenceFasta} .
            gunzip $REF_FILENAME
            REF_FILENAME="${REF_FILENAME%.gz}"
        else
            ln -s ~{referenceFasta}
        fi
        samtools faidx $REF_FILENAME

        # initialize script
        cmd=( /opt/dipcall/dipcall.kit/run-dipcall )

        # male samples need PAR region excluded
        if [[ ~{isMaleSample} == true ]]; then
            if [[ ~{referenceIsHS38} ]]; then
                cmd+=( -x /opt/dipcall/dipcall.kit/hs38.PAR.bed )
            elif [ -f "~{custom_PAR_bed}" ]; then
                cmd+=( -x "~{custom_PAR_bed}" )
            else
                echo "WARNING: Assuming reference is hg37 because input is not hg38 and no custom PAR provided."
                echo "Using built-in PAR regions for hg37 to exclude from calling..."
                cmd+=( -x /opt/dipcall/dipcall.kit/hs37d5.PAR.bed )
            fi

        fi

        # finalize script
        cmd+=( $PREFIX.dipcall/$PREFIX )
        cmd+=( $REF_FILENAME )
        cmd+=( $PAT_FILENAME )
        cmd+=( $MAT_FILENAME )

        # generate makefile
        "${cmd[@]}" >$PREFIX.mak

        # run dipcall
        make -j 2 -f $PREFIX.mak

        # finalize
        rm $PREFIX.dipcall/*sam.gz
        tar czvf $PREFIX.dipcall.tar.gz $PREFIX.dipcall/
        cp $PREFIX.dipcall/$PREFIX.dip.bed $PREFIX.dipcall.bed
        cp $PREFIX.dipcall/$PREFIX.dip.vcf.gz $PREFIX.dipcall.vcf.gz

        # cleanup
        rm $REF_FILENAME
        rm $MAT_FILENAME
        rm $PAT_FILENAME

	>>>
	output {
		File outputTarball = glob("*.dipcall.tar.gz")[0]
		File outputVCF = glob("*.dipcall.vcf.gz")[0]
		File outputBED = glob("*.dipcall.bed")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

# Task assumes contig names in PAR bed match the allosome names in Dipcall VCF, but otherwise name agnostic
# Also assumes that PAR are on the edges of X/Y with two connected components and nowhere else, as in humans
task CleanHaploidCalls {
    input {
        File dipcall_vcf
        File PAR_bed

        # Runtime arguments
        Int disk_size = 2*ceil(size(dipcall_vcf, "GB")) + 50
        Int cpu = 4
        Int memory = 8
    }

    command <<<
        set -xueo pipefail

        python << CODE

        import pysam
        from cyvcf2 import VCF, Writer


        input_file = "~{dipcall_vcf}"
        bed_file = "~{PAR_bed}"

        pysam.tabix_index(input_file, preset='vcf', force=True)

        vcf = VCF(input_file)
        vcf_out = Writer('cleaned.vcf.gz', vcf)

        # Function for formatting variants with missing genotypes in haploid regions
        def make_haploid_var(variant):
            gt_string = [x for x in variant.genotypes[0][:2] if (x != -1)]   # Grab non-missing entries; uses diploid assumption
            if len(gt_string) == 0:
                gt_string = '.'
            else:
                gt_string = gt_string[0]
            filter_string = variant.FILTER if variant.FILTER is not None else '.'
            var_string = '\t'.join([variant.CHROM, str(variant.POS), str(variant.ID), variant.REF,
                   ','.join(variant.ALT), f"{variant.QUAL:g}", filter_string, '.',
                   ':'.join(variant.FORMAT), f"{gt_string}:{','.join([str(x) for x in variant.format('AD').tolist()[0]])}"])

            return vcf_out.variant_from_string(var_string)

        # Process PAR bed file
        with open(bed_file) as bed:
            bed_text = bed.read().splitlines()

        regions = [x.split('\t') for x in bed_text]
        non_PAR = {}
        for i in range(0, len(regions), 2):
            non_PAR[regions[i][0]] = (regions[i][2], regions[i+1][1])

        # Iterate through variants and modify missing values in haploid regions
        for variant in vcf:
            # Check if variant in non-PAR X/Y region, i.e. expected haploid
            # Use > since bed coordinates are 0-based but VCF are 1-based
            if (variant.CHROM in set(non_PAR.keys())) and (variant.POS > int(non_PAR[variant.CHROM][0])) and (variant.POS <= int(non_PAR[variant.CHROM][1])):
                # Check if any missing genotypes
                if -1 in variant.genotypes[0]:
                    variant = make_haploid_var(variant)
            vcf_out.write_record(variant)

        vcf_out.close()

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/pysam:v1.1"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory + "GB"
    }

    output {
        File cleaned_vcf = "cleaned.vcf.gz"
    }
}

task IndexVCF {
    input {
        File input_vcf
        String sample_name = "dipcall"
    }

    command <<<
        set -xe

        echo "~{sample_name}" > samples.txt
        bcftools reheader -s samples.txt -o cleaned-indexed.vcf.gz ~{input_vcf}
        bcftools index -t -f cleaned-indexed.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.0"
        cpu: 2
        memory: "8GB"
        disks: "local-disk " + 100 + " HDD"
    }

    output {
        File vcf = "cleaned-indexed.vcf.gz"
        File vcf_index = "cleaned-indexed.vcf.gz.tbi"
    }
}
