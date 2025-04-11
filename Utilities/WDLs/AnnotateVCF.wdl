version 1.0

# A workflow for annotating VCF variant membership in an array of bed files
# Can also add GC content and sequence motif counts around each variant
workflow AnnotateVCF {
    input {
        File query_vcf
        File query_vcf_index
        File? truth_vcf
        File? truth_vcf_index
        File? truth_bed

        File ref_fasta
        File ref_fasta_index
        File fasta_dict

        Array[File] bed_files
        Array[String] bed_labels

        Boolean add_gc_content = true
        Int window_size = 25
        Boolean include_N_count = false
        Boolean annotate_sequence = false
        String? sequence_motif

        Array[String] extra_query_fields = []
        Array[String] extra_truth_fields = []

        Array[String] gatk_query_annotations = []
        Array[String] gatk_query_annotation_labels = []
        Array[String] gatk_truth_annotations = []
        Array[String] gatk_truth_annotation_labels = []
        File? query_bam
        File? query_bam_index
        File? gatk_jar
    }

    Array[String] gc_labels_1 = if (add_gc_content) then ["GC_CONTENT", "WINDOW_A", "WINDOW_C", "WINDOW_G", "WINDOW_T"] else []
    Array[String] gc_labels_2 = if (include_N_count && add_gc_content) then flatten([gc_labels_1, ["WINDOW_N"]]) else gc_labels_1
    Array[String] gc_labels_3 = if (annotate_sequence && add_gc_content) then flatten([gc_labels_2, ["SEQUENCE"]]) else gc_labels_2
    Array[String] gc_labels = if (defined(sequence_motif) && add_gc_content) then flatten([gc_labels_3, ["MOTIF_COUNTS"]]) else gc_labels_3

    # First annotate with GATK VariantAnnotator because it is most computationally expensive; for call caching
    if (length(gatk_truth_annotations) > 0 && defined(query_bam) && defined(query_bam_index) && defined(truth_vcf)) {
        call AnnotateVcfWithGatk as AnnotateTruthGatk {
            input:
                vcf=select_first([truth_vcf]),
                vcf_index=select_first([truth_vcf_index]),
                bam=select_first([query_bam]),
                bam_index=select_first([query_bam_index]),
                ref_fasta=ref_fasta,
                ref_fasta_index=ref_fasta_index,
                ref_fasta_dict=fasta_dict,
                annotations=gatk_truth_annotations,
                gatk_jar=gatk_jar,
                output_name="truth_gatk_annotated",
        }
    }

    if (length(gatk_query_annotations) > 0 && defined(query_bam) && defined(query_bam_index)) {
        call AnnotateVcfWithGatk as AnnotateQueryGatk {
            input:
                vcf=query_vcf,
                vcf_index=query_vcf_index,
                bam=select_first([query_bam]),
                bam_index=select_first([query_bam_index]),
                ref_fasta=ref_fasta,
                ref_fasta_index=ref_fasta_index,
                ref_fasta_dict=fasta_dict,
                annotations=gatk_query_annotations,
                gatk_jar=gatk_jar,
                output_name="query_gatk_annotated"
        }
    }

    # Next annotate using truth labels using vcfeval
    if (defined(truth_vcf)) {
        call VcfEvalAnnotate {
            input:
                query_vcf=select_first([AnnotateQueryGatk.annotated_vcf, query_vcf]),
                query_vcf_index=select_first([AnnotateQueryGatk.annotated_vcf_index, query_vcf_index]),
                truth_vcf=select_first([AnnotateTruthGatk.annotated_vcf, truth_vcf]),
                truth_vcf_index=select_first([AnnotateTruthGatk.annotated_vcf_index, truth_vcf_index]),
                truth_bed=truth_bed,
                ref_fasta=ref_fasta,
                ref_fasta_index=ref_fasta_index,
                bed_labels=bed_labels
        }

        call AnnotateVcfRegions as AnnotateTruthRegions {
            input:
                input_vcf=VcfEvalAnnotate.base_vcf,
                input_vcf_index=VcfEvalAnnotate.base_vcf_index,
                bed_files=bed_files,
                bed_labels=bed_labels
        }

        if (add_gc_content) {
            call AnnotateVcfGcContent as AnnotateTruthGcContent {
                input:
                    vcf=AnnotateTruthRegions.annotated_vcf,
                    vcf_index=AnnotateTruthRegions.annotated_vcf_index,
                    ref_fasta=ref_fasta,
                    ref_fasta_index=ref_fasta_index,
                    window_size=window_size,
                    include_N_count=include_N_count,
                    annotate_sequence=annotate_sequence,
                    sequence_motif=sequence_motif
            }
        }

        call ExtractData as ExtractTruthData {
            input:
                vcf=select_first([AnnotateTruthGcContent.annotated_vcf, AnnotateTruthRegions.annotated_vcf]),
                vcf_index=select_first([AnnotateTruthGcContent.annotated_vcf_index, AnnotateTruthRegions.annotated_vcf_index]),
                vcfeval_label="BASE",
                bed_labels=bed_labels,
                gc_labels=gc_labels,
                gatk_annotation_labels=gatk_truth_annotation_labels,
                extra_fields=extra_truth_fields
        }

    }

    call AnnotateVcfRegions as AnnotateQueryRegions {
        input:
            input_vcf=select_first([VcfEvalAnnotate.calls_vcf, AnnotateQueryGatk.annotated_vcf, query_vcf]),
            input_vcf_index=select_first([VcfEvalAnnotate.calls_vcf_index, AnnotateQueryGatk.annotated_vcf_index, query_vcf_index]),
            bed_files=bed_files,
            bed_labels=bed_labels
    }

    if (add_gc_content) {
        call AnnotateVcfGcContent as AnnotateQueryGcContent {
            input:
                vcf=AnnotateQueryRegions.annotated_vcf,
                vcf_index=AnnotateQueryRegions.annotated_vcf_index,
                ref_fasta=ref_fasta,
                ref_fasta_index=ref_fasta_index,
                window_size=window_size,
                include_N_count=include_N_count,
                annotate_sequence=annotate_sequence,
                sequence_motif=sequence_motif
        }
    }

    String query_label = if (defined(truth_vcf)) then "CALL" else "QUERY"

    call ExtractData as ExtractQueryData {
        input:
            vcf=select_first([AnnotateQueryGcContent.annotated_vcf, AnnotateQueryRegions.annotated_vcf]),
            vcf_index=select_first([AnnotateQueryGcContent.annotated_vcf_index, AnnotateQueryRegions.annotated_vcf_index]),
            vcfeval_label=query_label,
            bed_labels=bed_labels,
            gc_labels=gc_labels,
            gatk_annotation_labels=gatk_query_annotation_labels,
            extra_fields=extra_query_fields
    }


    output {
        File query_snp_table = ExtractQueryData.snp_table
        Int query_snp_count = ExtractQueryData.snp_count
        File query_indel_table = ExtractQueryData.indel_table
        Int query_indel_count = ExtractQueryData.indel_count

        File? truth_snp_table = ExtractTruthData.snp_table
        Int? truth_snp_count = ExtractTruthData.snp_count
        File? truth_indel_table = ExtractTruthData.indel_table
        Int? truth_indel_count = ExtractTruthData.indel_count
    }
}

task AnnotateVcfRegions {
    input {
        File input_vcf
        File input_vcf_index
        Array[File] bed_files
        Array[String] bed_labels

        Int disk_size = ceil(2 * size(input_vcf, "GB")) + 100
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xeuo pipefail

        bed_files=( ~{sep=" " bed_files} )
        bed_labels=( ~{sep=" " bed_labels} )

        # Move input VCF to avoid overwriting in loop
        mv ~{input_vcf} input.vcf.gz
        mv ~{input_vcf_index} input.vcf.gz.tbi

        # Annotate VCF with each BED file and corresponding label
        for i in ${!bed_files[@]}; do
            bed_file=${bed_files[$i]}
            bed_label=${bed_labels[$i]}

            bgzip $bed_file
            tabix -p bed "${bed_file}.gz"
            bcftools annotate input.vcf.gz \
                -a "${bed_file}.gz" \
                -c CHROM,POS \
                -h <(echo "##INFO=<ID=${bed_label},Number=0,Type=Flag,Description=\"Contained in ${bed_label}\">") \
                -m "+${bed_label}" \
                -o annotated.vcf.gz
            mv annotated.vcf.gz input.vcf.gz
            bcftools index -t input.vcf.gz
        done

        # Rename final output
        mv input.vcf.gz annotated_output.vcf.gz
        mv input.vcf.gz.tbi annotated_output.vcf.gz.tbi
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
        preemptible: 3
    }

    output {
        File annotated_vcf = "annotated_output.vcf.gz"
        File annotated_vcf_index = "annotated_output.vcf.gz.tbi"
    }
}


# Use bedtools nuc to annotate a GC window around VCF
# Use bedtools slop to create the windows around each variant
task AnnotateVcfGcContent {
    input {
        File vcf
        File vcf_index
        File ref_fasta
        File ref_fasta_index

        Int window_size = 25    # Bases to add to each side of variant
        Boolean include_N_count = false
        Boolean annotate_sequence = false
        String? sequence_motif    # Motif to count instances of in window
    }

    String annotate_sequence_command = if (annotate_sequence) then "-seq" else ""

    command <<<
        set -xeuo pipefail

        # Create a BED file with window stats around each variant
        bcftools query -f '%CHROM\t%POS\t%POS\n' ~{vcf} | \
            bedtools slop -i - -g ~{ref_fasta_index} -r ~{window_size} -l ~{window_size+1} | \
            bedtools nuc -fi ~{ref_fasta} -bed - ~{annotate_sequence_command} ~{"-C -pattern " + sequence_motif} > gc_content.bed

        # Annotate VCF with GC content
        bedtools slop -i gc_content.bed -g ~{ref_fasta_index} -b -~{window_size} > gc_content_restricted.bed
        bgzip gc_content_restricted.bed
        tabix -p bed gc_content_restricted.bed.gz

        echo '##INFO=<ID=GC_CONTENT,Number=1,Type=Float,Description="GC Content in window of size ~{window_size} around variant.">' > header.txt
        echo '##INFO=<ID=WINDOW_A,Number=1,Type=Integer,Description="Number of A bases in window of size ~{window_size} around variant.">' >> header.txt
        echo '##INFO=<ID=WINDOW_C,Number=1,Type=Integer,Description="Number of C bases in window of size ~{window_size} around variant.">' >> header.txt
        echo '##INFO=<ID=WINDOW_G,Number=1,Type=Integer,Description="Number of G bases in window of size ~{window_size} around variant.">' >> header.txt
        echo '##INFO=<ID=WINDOW_T,Number=1,Type=Integer,Description="Number of T bases in window of size ~{window_size} around variant.">' >> header.txt

        if ~{include_N_count}; then
            echo '##INFO=<ID=WINDOW_N,Number=1,Type=Integer,Description="Number of N bases in window of size ~{window_size} around variant.">' >> header.txt
            WINDOW_N=",WINDOW_N"
        else
            WINDOW_N=",-"
        fi

        if ~{annotate_sequence}; then
            echo '##INFO=<ID=SEQUENCE,Number=1,Type=String,Description="Sequence in window of size ~{window_size} around variant.">' >> header.txt
            SEQUENCE=",SEQUENCE"
        else
            SEQUENCE=""
        fi

        if ~{defined(sequence_motif)}; then
            echo '##INFO=<ID=MOTIF_COUNTS,Number=1,Type=Integer,Description="Number of instances of ~{sequence_motif} in window of size ~{window_size} around variant.">' >> header.txt
            MOTIF_COUNTS=",MOTIF_COUNTS"
        else
            MOTIF_COUNTS=""
        fi

        bcftools annotate ~{vcf} \
            -a gc_content_restricted.bed.gz \
            -c "CHROM,POS,-,-,GC_CONTENT,WINDOW_A,WINDOW_C,WINDOW_G,WINDOW_T${WINDOW_N},-,-${SEQUENCE}${MOTIF_COUNTS}" \
            -h header.txt \
            -o annotated.vcf.gz
        bcftools index -t annotated.vcf.gz

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
        disks: "local-disk 200 HDD"
        memory: "16 GB"
        cpu: 4
        preemptible: 3
    }

    output {
        File annotated_vcf = "annotated.vcf.gz"
        File annotated_vcf_index = "annotated.vcf.gz.tbi"
    }
}

# Use vcfeval to annotate input VCF with TP/FP labels using a truth VCF
task VcfEvalAnnotate {
    input {
        File query_vcf
        File query_vcf_index

        File truth_vcf
        File truth_vcf_index
        File? truth_bed

        File ref_fasta
        File ref_fasta_index

        Array[String] bed_labels
    }


    command <<<
        set -xueo pipefail

        rtg format -o rtg_ref ~{ref_fasta}
        rtg vcfeval \
            -b ~{truth_vcf} \
            -c ~{query_vcf} \
            -t rtg_ref \
            -o vcfeval_output \
            ~{"-e " + truth_bed} \
            --decompose \
            --output-mode annotate

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
        disks: "local-disk 300 HDD"
        memory: "16 GB"
        cpu: 8
    }

    output {
        File calls_vcf = "vcfeval_output/calls.vcf.gz"
        File calls_vcf_index = "vcfeval_output/calls.vcf.gz.tbi"
        File base_vcf = "vcfeval_output/baseline.vcf.gz"
        File base_vcf_index = "vcfeval_output/baseline.vcf.gz.tbi"
    }
}

# Take a list of GATK VariantAnnotations, a VCF, and a BAM and output annotated VCF using GATK VariantAnnotator
task AnnotateVcfWithGatk {
    input {
        File vcf
        File vcf_index
        File bam
        File bam_index
        File ref_fasta
        File ref_fasta_index
        File ref_fasta_dict

        File? gatk_jar

        Array[String] annotations = []
        String output_name = "annotated"
    }

    String gatk_name = if (defined(gatk_jar)) then "java -jar " + gatk_jar else "gatk"

    command <<<
        set -xeuo pipefail

        # Convert from cram to bam if necessary
        if [[ $(basename ~{bam} .cram) != $(basename ~{bam}) ]]; then
            samtools view -b -T ~{ref_fasta} ~{bam} -o input.bam
            samtools index input.bam
        else
            ln -s ~{bam} input.bam
            ln -s ~{bam_index} input.bam.bai
        fi

        # Update VCF sample name to match the RG sample
        RG_SAMPLE=$(samtools view -H input.bam | grep '^@RG' | sed 's/.*SM:\([^\t]*\).*/\1/' | uniq)
        VCF_SAMPLE=$(bcftools query -l ~{vcf})
        echo -e "${VCF_SAMPLE}\t${RG_SAMPLE}" > sample_map.txt
        echo -e "${RG_SAMPLE}\t${VCF_SAMPLE}" > sample_map_reverse.txt
        bcftools reheader ~{vcf} -s sample_map.txt -o input.vcf.gz
        bcftools index -t input.vcf.gz

        # Annotate VCF with GATK VariantAnnotator
        ~{gatk_name} VariantAnnotator \
            -R ~{ref_fasta} \
            -V input.vcf.gz \
            -I input.bam \
            -A ~{sep=" -A " annotations} \
            -O ~{output_name}-sample.vcf.gz

        # Update sample name back to original
        bcftools reheader ~{output_name}-sample.vcf.gz -s sample_map_reverse.txt -o ~{output_name}.vcf.gz
        bcftools index -t ~{output_name}.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.6.0.0"
        disks: "local-disk 400 HDD"
        memory: "16 GB"
        cpu: 8
    }

    output {
        File annotated_vcf = "~{output_name}.vcf.gz"
        File annotated_vcf_index = "~{output_name}.vcf.gz.tbi"
    }
}

# Extract annotated columns from input vcf using bcftools norm and query
task ExtractData {
    input {
        File vcf
        File vcf_index

        String vcfeval_label
        Array[String] bed_labels
        Array[String] gc_labels
        Array[String] gatk_annotation_labels = []
        Array[String] extra_fields = []
    }

    String gc_tab = if (length(gc_labels) > 0) then "\\t" else ""
    String gatk_tab = if (length(gatk_annotation_labels) > 0) then "\\t" else ""
    String extra_tab = if (length(extra_fields) > 0) then "\\t" else ""

    String vcfeval_header_label = if (vcfeval_label == "QUERY") then "" else "~{vcfeval_label}\t"
    String bcftools_query_label = if (vcfeval_label == "QUERY") then "" else "~{vcfeval_header_label}%"
    String bcftools_exclude = if (vcfeval_label == "QUERY") then "" else "-e '~{vcfeval_label}=\"IGN\" || ~{vcfeval_label}=\"OUT\"'"

    command <<<
        set -xueo pipefail

        echo -e "CHROM\tPOS\tREF\tALT\t~{vcfeval_header_label}~{sep="\\t" bed_labels}~{gc_tab}~{sep="\\t" gc_labels}~{gatk_tab}~{sep="\\t" gatk_annotation_labels}~{extra_tab}~{sep="\\t" extra_fields}" > ~{vcfeval_label}_snp_table.tsv
        echo -e "CHROM\tPOS\tREF\tALT\t~{vcfeval_header_label}~{sep="\\t" bed_labels}~{gc_tab}~{sep="\\t" gc_labels}~{gatk_tab}~{sep="\\t" gatk_annotation_labels}~{extra_tab}~{sep="\\t" extra_fields}" > ~{vcfeval_label}_indel_table.tsv

        # Prepare extract strings
        GC_LABELS='~{true="\\t%" false="" length(gc_labels) > 0}~{sep="\\t%" gc_labels}'
        GATK_LABELS='~{true="\\t%" false="" length(gatk_annotation_labels) > 0}~{sep="\\t%" gatk_annotation_labels}'
        EXTRA_FIELDS='~{true="\\t%" false="" length(extra_fields) > 0}~{sep="\\t%" extra_fields}'
        QUERY_STRING="%CHROM\t%POS\t%REF\t%ALT\t%~{bcftools_query_label}~{sep="\\t%" bed_labels}${GC_LABELS}${GATK_LABELS}${EXTRA_FIELDS}\n"

        bcftools norm -m -any ~{vcf} -o normalized.vcf.gz
        bcftools index -t normalized.vcf.gz

        # Split into snp and indel tables
        bcftools view -v snps normalized.vcf.gz | bcftools query ~{bcftools_exclude} -f $QUERY_STRING >> ~{vcfeval_label}_snp_table.tsv
        bcftools view -v indels normalized.vcf.gz | bcftools query ~{bcftools_exclude} -f $QUERY_STRING >> ~{vcfeval_label}_indel_table.tsv

        # Compress tables
        gzip ~{vcfeval_label}_snp_table.tsv
        gzip ~{vcfeval_label}_indel_table.tsv

        # Count entries in each table
        zcat ~{vcfeval_label}_snp_table.tsv.gz | wc -l > snp_count.txt
        zcat ~{vcfeval_label}_indel_table.tsv.gz | wc -l > indel_count.txt
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
        disks: "local-disk 200 HDD"
        memory: "8 GB"
        cpu: 4
        preemptible: 3
    }

    output {
        File snp_table = "~{vcfeval_label}_snp_table.tsv.gz"
        Int snp_count = read_int("snp_count.txt")
        File indel_table = "~{vcfeval_label}_indel_table.tsv.gz"
        Int indel_count = read_int("indel_count.txt")
    }
}
