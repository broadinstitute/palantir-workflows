version 1.0

task cram_to_bam {
    input {
        File reference_fasta
        File reference_fai
        File reference_dict
        File cram_file
        String sample_id

        Int? cpu = 4
        Int? mem_gb = 32
        Int? disk_size_gb = 512
    }

    command <<<
        set -e
        set -o pipefail

        ln -vs ~{reference_fasta} reference.fasta
        ln -vs ~{reference_fai} reference.fasta.fai
        ln -vs ~{reference_dict} reference.dict

        samtools view -h -T reference.fasta ~{cram_file} |
        samtools view -b -o ~{sample_id}.bam -
        samtools index -b ~{sample_id}.bam
        mv ~{sample_id}.bam.bai ~{sample_id}.bai
    >>>

    output {
        File output_bam = "~{sample_id}.bam"
        File output_bam_index = "~{sample_id}.bai"
    }

    runtime {
        cpu: cpu
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
    }
}

task read_counter {
    input {
        File bam_file
        File bam_index
        String chrs
        String sample_id
        Int bin_size
        Int? min_qual = 20

        Int? cpu = 4
        Int? mem_gb = 32
        Int? disk_size_gb = ceil(size(bam_file, "GiB") + 100)
    }

    command <<<
        ln -vs ~{bam_index} ~{sample_id}.bam.bai
        ln -vs ~{bam_file} ~{sample_id}.bam

        /HMMcopy/bin/readCounter \
        ~{sample_id}.bam \
        -c ~{chrs} \
        -w ~{bin_size} \
        -q ~{min_qual} > ~{sample_id}.bin~{bin_size}.wig
    >>>

    runtime {
        cpu: cpu
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        docker: "us.gcr.io/tag-team-160914/bloodbiopsy-hmmcopy:0.0.1"
    }

    output {
        File wig_file = "~{sample_id}.bin~{bin_size}.wig"
    }
}

task ichor_cna_task {
    input {
        File wig_file
        File? normal_panel
        String genome_build
        Int bin_size_kb
        Float? mean_depth
        String sample_id
        String ploidy
        String normal
        Int max_cn
        Boolean include_homd
        String chrs
        String chr_train
        String chr_normalize
        String genome_style
        Boolean estimate_normal
        Boolean estimate_ploidy
        Boolean estimate_clonality
        String sc_states
        File centromere
        File? exons
        Float txn_e
        Int txn_strength
        Int min_segment_bins
        Float min_map_score
        Float frac_reads_chrY_male
        Float max_frac_cna_subclone
        Float max_frac_genome_subclone
        Float alt_frac_threshold
        Int lambda_scale_hyper_param
        Int rm_centromere_flank_length
        String? plot_file_type = "pdf"
        String plot_ylim
        File input_gc_wig
        File input_map_wig

        Int? cpu = 4
        Int? mem_gb = 32
        Int? disk_size_gb = 256
    }

    command <<<
        Rscript /runIchorCNA.R \
        --id ~{sample_id} \
        --outDir ./ --libdir /ichorCNA \
        --WIG ~{wig_file} \
        --gcWig  ~{input_gc_wig} \
        --mapWig ~{input_map_wig} \
        --normalPanel ~{default="None" normal_panel} \
        --ploidy "~{ploidy}" \
        --normal "~{normal}" \
        ~{"--coverage " + mean_depth} \
        --maxCN ~{max_cn} \
        --includeHOMD ~{include_homd} \
        --chrs "~{chrs}" \
        --chrTrain "~{chr_train}" \
        --chrNormalize "~{chr_normalize}" \
        --genomeStyle "~{genome_style}" \
        --genomeBuild "~{genome_build}" \
        --estimateNormal ~{estimate_normal} \
        --estimatePloidy ~{estimate_ploidy}  \
        --estimateScPrevalence ~{estimate_clonality} \
        --scStates "~{sc_states}" \
        --centromere ~{centromere} \
        ~{"--exons.bed " + exons} \
        --txnE ~{txn_e} \
        --txnStrength ~{txn_strength} \
        --minSegmentBins ~{min_segment_bins} \
        --minMapScore ~{min_map_score} \
        --lambdaScaleHyperParam ~{lambda_scale_hyper_param} \
        --fracReadsInChrYForMale ~{frac_reads_chrY_male} \
        --maxFracGenomeSubclone ~{max_frac_genome_subclone} \
        --altFracThreshold ~{alt_frac_threshold} \
        --maxFracCNASubclone ~{max_frac_cna_subclone} \
        --rmCentromereFlankLength ~{rm_centromere_flank_length} \
        --plotFileType ~{plot_file_type} \
        --plotYLim "~{plot_ylim}"

        # Zip optimal solutions
        mkdir ~{sample_id}.optimalSolution
        cp ~{sample_id}/~{sample_id}_genomeWide.pdf ~{sample_id}.cna.seg ~{sample_id}.seg.txt ~{sample_id}.seg ~{sample_id}.optimalSolution/
        zip -r ~{sample_id}.optimalSolution.zip ~{sample_id}.optimalSolution

        # Generate list of out solutions
        Rscript /gatherOutSolutions.R --id ~{sample_id} --rdata ~{sample_id}.RData
    >>>

    runtime {
        cpu: cpu
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        docker: "us.gcr.io/tag-team-160914/bloodbiopsy-ichorcna:0.2.1"
    }

    output {
        File corr_depth = "~{sample_id}.correctedDepth.txt"
        File params = "~{sample_id}.params.txt"
        File cna = "~{sample_id}.cna.seg"
        File seg_txt = "~{sample_id}.seg.txt"
        File seg = "~{sample_id}.seg"
        File rdata = "~{sample_id}.RData"

        File all_genomewide_plots = "~{sample_id}/~{sample_id}_genomeWide_all_sols.pdf"
        File bias = "~{sample_id}/~{sample_id}_bias.pdf"
        File tpdf = "~{sample_id}/~{sample_id}_tpdf.pdf"
        File correct = "~{sample_id}/~{sample_id}_correct.pdf"
        Array[File] per_chromosome_plots = glob("~{sample_id}/~{sample_id}_CNA*")

        File optimal_solution = "~{sample_id}.optimalSolution.zip"
        File out_solutions = "~{sample_id}.outSolutions.zip"
    }
}

task extract_ichor_params {
    input {
        File params
        Float max_frac_cna_subclone
        Float max_frac_genome_subclone

        Int? cpu = 1
        Int? mem_gb = 8
        Int? disk_size_gb = 64
    }

    command <<<
        cut -f1,2 ~{params} | grep -v ^$ | tr "\t" " " > params_table.txt

        python << CODE

        solutions = open('~{params}').readlines()
        sols = [x.rstrip("\n").split("\t") for x in solutions if x.startswith("n0.")]
        log_lik=0

        for sol in sols:
            if int(float(sol[6]))>log_lik and float(sol[4])<~{max_frac_genome_subclone} and float(sol[5])<~{max_frac_cna_subclone}:
                log_lik=int(float(sol[6]))

        params = open("params_table.txt","r").readlines()
        params = [x.rstrip("\n").split(": ") for x in params if ":" in x]
        param_dict = dict()

        for a,b in params:
            param_dict[a] = b

        p=open("tumor_fraction","w"); p.write(param_dict["Tumor Fraction"])
        p=open("ploidy","w"); p.write(param_dict["Ploidy"])
        p=open("subclone_fraction","w"); p.write(param_dict["Subclone Fraction"])
        p=open("fraction_genome_subclonal","w"); p.write(param_dict["Fraction Genome Subclonal"])
        p=open("fraction_cna_subclonal","w"); p.write(param_dict["Fraction CNA Subclonal"])
        p=open("gc-map_correction_mad","w"); p.write(param_dict["GC-Map correction MAD"])
        p=open("top_solution_log_likelihood","w"); p.write(str(log_lik))
        p.close()

        CODE
    >>>

    runtime {
        cpu: cpu
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        docker: "python:3"
    }

    output {
        Float tumor_fraction = read_float("tumor_fraction")
        Float ploidy = read_float("ploidy")
        String subclone_fraction = read_string("subclone_fraction")
        String fraction_genome_subclonal = read_string("fraction_genome_subclonal")
        String fraction_cna_subclonal = read_string("fraction_cna_subclonal")
        String gc_map_correction_mad = read_string("gc-map_correction_mad")
        Int top_solution_log_likelihood = read_int("top_solution_log_likelihood")
    }
}

task bundle_per_chromosome_plots {
    input {
        Array[File] chrom_plots
        String sample_id

        Int? cpu = 1
        Int? mem_gb = 8
        Int? disk_size_gb = 64
    }

    command <<<
        set -e
        CHROM_PLOTS=`ls ~{sep=" " chrom_plots} | sort -V`

        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite \
        -sOutputFile="~{sample_id}_OptimalSolutionPerChrom.pdf" \
        $CHROM_PLOTS
    >>>

    runtime {
        cpu: cpu
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        docker: "us.gcr.io/tag-team-160914/tag-tools:1.0.0"
    }

    output {
        File output_plot = "~{sample_id}_OptimalSolutionPerChrom.pdf"
    }
}

workflow ichor_cna_cram {
    input {
        File cram_file
        File cram_index
        File reference_fasta
        File reference_fai
        File reference_dict
        String sample_id
        String genome_build
        File? normal_panel
        Int bin_size_kb
        Int bin_size = bin_size_kb * 1000
        File centromere
        File gc_wig
        File map_wig
        String genome_style
        String chr_counter      # All chromosomes for counting reads, as a list
        String chrs             # Autosomal and female chromosome for running ichor, as an R command
        String chr_train        # Autosomal chromosomes to estimate ichor params, as an R command
        Float max_frac_cna_subclone
        Float max_frac_genome_subclone
    }

    call cram_to_bam {
        input:
            reference_fasta = reference_fasta,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            cram_file = cram_file,
            sample_id = sample_id
    }

    call read_counter {
        input:
            bam_file = cram_to_bam.output_bam,
            bam_index = cram_to_bam.output_bam_index,
            sample_id = sample_id,
            bin_size = bin_size,
            chrs = chr_counter
    }

    call ichor_cna_task {
        input:
            wig_file = read_counter.wig_file,
            genome_build = genome_build,
            sample_id = sample_id,
            normal_panel = normal_panel,
            bin_size_kb = bin_size_kb,
            chrs = chrs,
            chr_train = chr_train,
            centromere = centromere,
            genome_style = genome_style,
            input_gc_wig = gc_wig,
            input_map_wig = map_wig,
            max_frac_cna_subclone = max_frac_cna_subclone,
            max_frac_genome_subclone = max_frac_genome_subclone
    }

    call extract_ichor_params {
        input:
            params = ichor_cna_task.params,
            max_frac_cna_subclone = max_frac_cna_subclone,
            max_frac_genome_subclone = max_frac_genome_subclone
    }

    call bundle_per_chromosome_plots {
        input:
            chrom_plots = ichor_cna_task.per_chromosome_plots,
            sample_id = sample_id
    }

    output {
        File wig_file = read_counter.wig_file

        File all_genomewide_plots = ichor_cna_task.all_genomewide_plots
        File corr_depth = ichor_cna_task.corr_depth
        File cna = ichor_cna_task.cna
        File seg_txt = ichor_cna_task.seg_txt
        File seg = ichor_cna_task.seg
        File rdata = ichor_cna_task.rdata

        Float tumor_fraction = extract_ichor_params.tumor_fraction
        Float ploidy = extract_ichor_params.ploidy
        String subclone_fraction = extract_ichor_params.subclone_fraction
        String fraction_genome_subclonal = extract_ichor_params.fraction_genome_subclonal
        String fraction_cna_subclonal = extract_ichor_params.fraction_cna_subclonal
        String gc_map_correction_mad = extract_ichor_params.gc_map_correction_mad
        Int top_solution_log_likelihood = extract_ichor_params.top_solution_log_likelihood

        File bias = ichor_cna_task.bias
        File tpdf = ichor_cna_task.tpdf
        File correct = ichor_cna_task.correct
        File params = ichor_cna_task.params

        File optimal_solution = ichor_cna_task.optimal_solution
        File out_solutions = ichor_cna_task.out_solutions
        File per_chromosome_plots = bundle_per_chromosome_plots.output_plot
    }
}