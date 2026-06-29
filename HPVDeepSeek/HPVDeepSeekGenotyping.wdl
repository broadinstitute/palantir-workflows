version 1.0

task FastQC {
    input {
        File r1_fastq
        File r2_fastq
        Int fastqc_thread_memory = 4096

        Int cpu = 2
        Int num_threads = 2
        Int memory_gb = 16
        Int disk_size_gb = ceil((2 * (size(r1_fastq, "GiB") + size(r2_fastq, "GiB"))) + 50)
        Int min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    String r1_fastq_name = sub(sub(sub(sub(basename(r1_fastq), "\\.fastq.gz$", ""), "\\.fq.gz$", ""), "\\.fastq$", ""), "\\.fq$", "")
    String r2_fastq_name = sub(sub(sub(sub(basename(r2_fastq), "\\.fastq.gz$", ""), "\\.fq.gz$", ""), "\\.fastq$", ""), "\\.fq$", "")

    command <<<
        fastqc ~{r1_fastq} ~{r2_fastq} \
        --threads ~{num_threads} \
        --memory ~{fastqc_thread_memory} \
        --outdir .

        unzip -p ~{r1_fastq_name}_fastqc.zip ~{r1_fastq_name}_fastqc/fastqc_data.txt | gzip > ~{r1_fastq_name}.fastqc_data.txt.gz
        unzip -p ~{r2_fastq_name}_fastqc.zip ~{r2_fastq_name}_fastqc/fastqc_data.txt | gzip > ~{r2_fastq_name}.fastqc_data.txt.gz
    >>>

    output {
        File r1_fastqc_html = "~{r1_fastq_name}_fastqc.html"
        File r1_fastqc_zip =  "~{r1_fastq_name}_fastqc.zip"
        File r1_fastqc_data = "~{r1_fastq_name}.fastqc_data.txt.gz"
        File r2_fastqc_html = "~{r2_fastq_name}_fastqc.html"
        File r2_fastqc_zip =  "~{r2_fastq_name}_fastqc.zip"
        File r2_fastqc_data = "~{r2_fastq_name}.fastqc_data.txt.gz"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
    }
}

task FastqToUbam {
    input {
        String output_basename
        File r1_fastq
        File r2_fastq
        String read_group_id
        String read_group_sample_name
        String read_group_library_name
        String read_group_platform
        String read_group_platform_unit
        String read_group_description

        Int cpu = 1
        Int memory_gb = 16
        Int disk_size_gb = ceil((2.5 * (size(r1_fastq, "GiB") + size(r2_fastq, "GiB"))) + 50)
        Int min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    command <<<
        gatk FastqToSam \
        --FASTQ ~{r1_fastq} \
        --FASTQ2 ~{r2_fastq} \
        --OUTPUT ~{output_basename}.unmapped.bam \
        --READ_GROUP_NAME ~{read_group_id} \
        --SAMPLE_NAME ~{read_group_sample_name} \
        --LIBRARY_NAME ~{read_group_library_name} \
        --PLATFORM ~{read_group_platform} \
        --PLATFORM_UNIT ~{read_group_platform_unit} \
        --DESCRIPTION ~{read_group_description}
    >>>

    output {
        File ubam = "~{output_basename}.unmapped.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# Extract the UMI sequence from the first 3 bases of each read, skip the next 2 bases
# and add the resulting UMI to the RX tag and read name in the output BAM.
#
# --molecular-index-tags RX: Tells fgbio to place the extracted UMI into the RX tag of each read.
# --read-structure 3M2S+T 3M2S+T: Specifies the regular expressions used to extract the UMI from each read sequence:
# 3M2S+T means:
# 3M: Match the first 3 bases (UMI).
# 2S: Skip (soft clip) the next 2 bases.
# +T: The rest of the read is the template
# This regex is applied to both R1 and R2 (paired-end reads).
# --annotate-read-names true: Indicates that the UMI should be appended to the read name (QNAME) in addition to the RX tag.
task ExtractUMIs {
    input {
        String output_basename
        File input_ubam
        String read_group_tag = "RX"
        String read_structure
        String append_umi_to_qname = "true"

        Int cpu = 1
        Int memory_gb = 16
        Int disk_size_gb = ceil((2.5 * size(input_ubam, "GiB")) + 50)
        Int min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    command <<<
        fgbio ExtractUmisFromBam \
        --input ~{input_ubam} \
        --output ~{output_basename}.umi_extracted.unmapped.bam \
        --read-structure ~{read_structure} ~{read_structure} \
        --molecular-index-tags ~{read_group_tag} \
        --annotate-read-names ~{append_umi_to_qname}
    >>>

    output {
        File umi_extracted_bam = "~{output_basename}.umi_extracted.unmapped.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# Convert UMI extracted uBAM to FASTQ before adapter trimming and filtering
task UmiExtractedBamToFastq {
    input {
        File umi_extracted_bam
        String output_basename

        Int cpu = 1
        Int memory_gb = 16
        Int disk_size_gb = ceil((2.5 * size(umi_extracted_bam, "GiB")) + 50)
        Int min_ssd_size_gb = 1024
        Boolean use_ssd = true
    }

    command <<<
        gatk SamToFastq \
        --INPUT ~{umi_extracted_bam} \
        --FASTQ ~{output_basename}_R1.fastq \
        --SECOND_END_FASTQ ~{output_basename}_R2.fastq
    >>>

    output {
        File umi_extracted_fastq1 = "~{output_basename}_R1.fastq"
        File umi_extracted_fastq2 = "~{output_basename}_R2.fastq"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# Clean and filter paired-end reads for downstream analysis,
# improving quality by trimming adapters, low-quality bases, and short reads.
# Filtering and Trimming Options
# -g: Detect and remove adapter sequences automatically.
# -W 5: Sliding window size of 5 for quality filtering.
# -q 20: Minimum quality score to keep during sliding window.
# -u 40: Remove reads with >40% low-quality bases.
# -3: Trim poly-G tails (common in Illumina NextSeq/NovaSeq).
# -l 75: Discard reads shorter than 75 bases after trimming.
# -c: Trim adapters even if only one end matches.
task TrimAndFilter {
    input {
        File fastq1
        File fastq2
        String output_basename

        Int cpu = 3
        Int memory_gb = 16
        Int disk_size_gb = ceil((2.5 * (size(fastq1, "GiB") + size(fastq2, "GiB"))) + 50)
        Int min_ssd_size_gb = 1024
        Boolean use_ssd = true
    }

    command <<<
        fastp \
        -i ~{fastq1} \
        -I ~{fastq2} \
        -o ~{output_basename}_R1.trimmed.fastq \
        -O ~{output_basename}_R2.trimmed.fastq \
        -g \
        -W 5 \
        -q 20 \
        -u 40 \
        -3 \
        -l 75 \
        -c \
        -h ~{output_basename}.fastp_report.html \
        -j ~{output_basename}.fastp_report.json
    >>>

    output {
        File fastq1_trimmed = "~{output_basename}_R1.trimmed.fastq"
        File fastq2_trimmed = "~{output_basename}_R2.trimmed.fastq"
        File fastp_report_html = "~{output_basename}.fastp_report.html"
        File fastp_report_json = "~{output_basename}.fastp_report.json"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# Align paired-end sequencing reads to a reference genome using BWA-MEM
# -R '@RG\tID:A\tDS:KAPA_TE\tPL:ILLUMINA\tLB:lib1\tSM:sample1\tPU:unit1'\-M {reference_genome} \
# -t {thread_count}: Number of threads to use for parallel processing.
# -K 100000000: Increases the internal buffer size for performance (100M input bytes per batch).
# -R '@RG\tID:A\tDS:KAPA_TE\tPL:ILLUMINA\tLB:lib1\tSM:sample1\tPU:unit1': Specifies the read group info in SAM header:
# ID:A: Identifier for this read group.
# DS:KAPA_TE: Description, here a prep kit name.
# PL:ILLUMINA: Platform (e.g., Illumina).
# LB:lib1: Library ID.
# SM:sample1: Sample name.
# PU:unit1: Platform unit, often includes lane/barcode info.
# -M: Marks shorter split hits as secondary (needed for compatibility with Picard and GATK).
task BwaMem {
    input {
        String output_basename
        File fastq1
        File fastq2
        File reference
        File bwa_idx_amb # !UnusedDeclaration
        File bwa_idx_ann # !UnusedDeclaration
        File bwa_idx_bwt # !UnusedDeclaration
        File bwa_idx_pac # !UnusedDeclaration
        File bwa_idx_sa  # !UnusedDeclaration
        String read_group_id
        String read_group_sample_name
        String read_group_library_name
        String read_group_platform
        String read_group_platform_unit
        String read_group_description
        Boolean soft_clip_supplementary_alignments = false

        Int cpu = 32
        Int num_threads = 32
        Int memory_gb = 64
        Int disk_size_gb = ceil((3 * (size(fastq1, "GiB") + size(fastq2, "GiB"))) + size(reference, "GiB") + 100)
        Int min_ssd_size_gb = 1024
        Boolean use_ssd = true
    }

    String supplementary_alignment_clipping_option = if soft_clip_supplementary_alignments then "-Y" else ""

    command <<<
        bwa mem \
        -t ~{num_threads} \
        -K 100000000 \
        -R '@RG\tID:~{read_group_id}\tDS:~{read_group_description}\tPL:~{read_group_platform}\tLB:~{read_group_library_name}\tSM:~{read_group_sample_name}\tPU:~{read_group_platform_unit}' \
        ~{supplementary_alignment_clipping_option} \
        -M \
        ~{reference} \
        ~{fastq1} \
        ~{fastq2} \
        | samtools view --threads ~{num_threads} -o ~{output_basename}.bam -
    >>>

    output {
        File bam = "~{output_basename}.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# Sort and index aligned BAM
task SortAndIndexBam {
    input {
        File bam
        String samtools_thread_memory = "1024M"

        Int cpu = 7
        Int num_threads = 7
        Int memory_gb = 16
        Int disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
        Int min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    String prefix = basename(bam, ".bam")

    command <<<
        samtools sort \
        -o ~{prefix}.sorted.bam \
        -O bam \
        -T ~{prefix}.bam.temp \
        -@ ~{num_threads} \
        -m ~{samtools_thread_memory} \
        ~{bam}

        samtools index ~{prefix}.sorted.bam
    >>>

    output {
        File sorted_bam = "~{prefix}.sorted.bam"
        File sorted_bam_index = "~{prefix}.sorted.bam.bai"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# GATK Sort BAM by queryname
task GATKSortBam {
    input {
        File bam

        Int cpu = 1
        Int memory_gb = 16
        Int disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
        Int min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    String prefix = basename(bam, ".bam")

    command <<<
        gatk SortSam \
        --INPUT ~{bam} \
        --OUTPUT ~{prefix}.sorted.bam \
        --SORT_ORDER queryname
    >>>

    output {
        File sorted_bam = "~{prefix}.sorted.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

task MergeBamAlignment {
    input {
        String output_basename
        File aligned_bam
        File unmapped_bam
        File reference
        File reference_fai
        File reference_dict
        String? extra_args

        Int cpu = 1
        Int memory_gb = 16
        Int disk_size_gb = ceil((2.5 * size(aligned_bam, "GiB") + size(unmapped_bam, "GiB")) + 100)
        Int min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    command <<<
        gatk MergeBamAlignment \
        --ALIGNED_BAM ~{aligned_bam} \
        --UNMAPPED_BAM ~{unmapped_bam} \
        --OUTPUT ~{output_basename}.bam \
        --REFERENCE_SEQUENCE ~{reference} \
        ~{extra_args}
    >>>

    output {
        File merged_bam = "~{output_basename}.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

task FilterAndGroupReadsByUMI {
    input {
        String output_basename
        File merged_bam
        File reference # !UnusedDeclaration
        File reference_fai
        File reference_dict
        Boolean is_duplex

        Int cpu = 1
        Int memory_gb = 16
        Int disk_size_gb = ceil((2.5 * size(merged_bam, "GiB")) + 100)
        Int min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    String strategy = if is_duplex then "paired" else "adjacency"
    String output_type = if is_duplex then "duplex" else "simplex"

    command <<<
        samtools view -f 2 -q 1 -bh ~{merged_bam} -o ~{output_basename}.merged.filtered.bam

        fgbio GroupReadsByUmi \
        --input ~{output_basename}.merged.filtered.bam \
        --output ~{output_basename}.~{output_type}.umi_grouped.bam \
        --strategy ~{strategy} \
        --edits 1 \
        --raw-tag RX \
        --family-size-histogram ~{output_basename}.~{output_type}.umi_group_data.txt
    >>>

    output {
        File umi_grouped_bam = "~{output_basename}.~{output_type}.umi_grouped.bam"
        File umi_group_data = "~{output_basename}.~{output_type}.umi_group_data.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# Collapse PCR duplicates with the same UMI and alignment coordinates into a single consensus read, improving accuracy by reducing sequencing errors.
# --input: BAM file with UMI-grouped reads (contains MI tags).
# --output: BAM with consensus reads, still unmapped (requires re-alignment).
# --error-rate-post-umi 40: Error rate assumed after UMI attachment (Phred scale = 1 in 10,000).
# --error-rate-pre-umi 45: Error rate before UMI, typically higher (1 in 32,000).
# --output-per-base-tags false: Disables additional tags with per-base information (smaller output).
# --min-reads 1: Minimum number of reads per group to form a consensus (even singletons retained).
# --max-reads 50: Caps number of reads used per group (avoids long runtimes).
# --min-input-base-quality 20: Filters out low-quality bases before building consensus.
# --read-name-prefix='consensus': Prefix for read names in output BAM.
task CallMolecularConsensusReads {
    input {
        String output_basename
        File umi_grouped_bam
        String read_group_id

        Int cpu = 1
        Int memory_gb = 16
        Int disk_size_gb = ceil((2.5 * size(umi_grouped_bam, "GiB")) + 100)
    }

    command <<<
        fgbio CallMolecularConsensusReads \
        --input ~{umi_grouped_bam} \
        --output ~{output_basename}.simplex.umi_consensus.unmapped.bam \
        --error-rate-post-umi 40 \
        --error-rate-pre-umi 45 \
        --output-per-base-tags false \
        --min-reads 1 \
        --max-reads 50 \
        --min-input-base-quality 20 \
        --read-name-prefix 'consensus' \
        --read-group-id ~{read_group_id}
    >>>

    output {
        File umi_consensus_unmapped_bam = "~{output_basename}.simplex.umi_consensus.unmapped.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

task CallDuplexConsensusReads {
    input {
        String output_basename
        File umi_grouped_bam
        String read_group_id

        Int cpu = 1
        Int memory_gb = 16
        Int disk_size_gb = ceil((2.5 * size(umi_grouped_bam, "GiB")) + 100)
    }

    command <<<
        fgbio CallDuplexConsensusReads \
        --input ~{umi_grouped_bam} \
        --output ~{output_basename}.duplex.umi_consensus.unmapped.bam \
        --error-rate-post-umi 40 \
        --error-rate-pre-umi 45 \
        --min-reads 1 \
        --max-reads-per-strand 50 \
        --min-input-base-quality 20 \
        --read-name-prefix 'consensus' \
        --read-group-id ~{read_group_id} \
        --consensus-call-overlapping-bases true \
        --stats ~{output_basename}.duplex.consensus.stats.txt
    >>>

    output {
        File umi_consensus_unmapped_bam = "~{output_basename}.duplex.umi_consensus.unmapped.bam"
        File consensus_stats = "~{output_basename}.duplex.consensus.stats.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/fgbio@sha256:b6869a0ae243d9f1b183e4a986fbe0853df2a56a1c6d7c0fec2965b6d8a7af1d"
    }
}

# Convert consensus BAM to FASTQ:
task ConsensusBamToFastq {
    input {
        String output_basename
        File umi_consensus_unmapped_bam

        Int cpu = 2
        Int memory_gb = 16
        Int disk_size_gb = ceil((2.5 * size(umi_consensus_unmapped_bam, "GiB")) + 50)
    }

    command <<<
        gatk SamToFastq \
        --INPUT ~{umi_consensus_unmapped_bam} \
        --FASTQ ~{output_basename}_R1.consensus.fastq \
        --SECOND_END_FASTQ ~{output_basename}_R2.consensus.fastq
    >>>

    output {
        File consensus_unmapped_fastq1 = "~{output_basename}_R1.consensus.fastq"
        File consensus_unmapped_fastq2 = "~{output_basename}_R2.consensus.fastq"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# Reads and Coverage Calculation, and results generation
task SamtoolsCoverage {
    input {
        String output_basename
        File bam

        Int cpu = 2
        Int memory_gb = 16
        Int disk_size_gb = ceil((2 * size(bam, "GiB")) + 50)
    }

    command <<<
        samtools coverage ~{bam} > ~{output_basename}.coverage.txt
    >>>

    output {
        File coverage = "~{output_basename}.coverage.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# HPV+ Classification
# Sample is considered HPV+ if the following thresholds are met:
# Duplex read count ≥ 2
task DetermineHPVStatus {
    input {
        String output_basename
        File coverage
        File low_risk_hpv_genotypes

        Int cpu = 1
        Int memory_gb = 8
        Int disk_size_gb = 32
    }

    command <<<
        set -e
        python3 <<CODE

        import pandas as pd

        low_risk_hpv_genotype_list = []
        with open("~{low_risk_hpv_genotypes}", 'r') as f:
            low_risk_hpv_genotype_list = f.read().splitlines()

        df = pd.read_csv("~{coverage}", sep = '\t')
        df = df.rename(columns = {"#rname": "rname"})
        df = df[["rname", "numreads", "coverage"]]

        df = df[(df["rname"].str.startswith("HPV")) & (df["numreads"] >= 2)]
        df["Is_Reportable"] = ~df.rname.isin(low_risk_hpv_genotype_list)

        df = df.rename(columns = {"rname": "HPV_Genotype", "numreads": "Num_Duplex_Reads", "coverage": "%_Genomic_Coverage"})
        df.to_csv("~{output_basename}.hpv_status.tsv", sep = '\t', index = False)

        CODE
    >>>

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/simple_pysam@sha256:f7f71cf1996056c32a4c7ad5ef6a855093383ab5c919a3539874e56a85539256"
    }

    output {
        File hpv_status = "~{output_basename}.hpv_status.tsv"
    }
}

# Human SNP Genotyping
task GenotypeSNPsHuman {
    input {
        String output_basename
        File bam
        File bai
        File human_snp_targets_bed
        File reference

        Int cpu = 2
        Int memory_gb = 16
        Int disk_size_gb = ceil((2.5 * size(bam, "GiB")) + 50)
    }

    command <<<
        bcftools mpileup -f ~{reference} -R ~{human_snp_targets_bed} ~{bam} 2> ~{output_basename}.mpileup.log \
        | bcftools call -mv -Ov -o ~{output_basename}.vcf 2> ~{output_basename}.call.log
    >>>

    output {
        File vcf = "~{output_basename}.vcf"
        File mpileup_log = "~{output_basename}.mpileup.log"
        File call_log = "~{output_basename}.call.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.4"
    }
}

task CollectAlignmentSummaryMetrics {
    input {
        File bam
        File bai
        File reference
        File reference_fai
        File reference_dict

        Int cpu = 2
        Int memory_gb = 16
        Int disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
    }

    String prefix = basename(bam, ".sorted.bam")

    command <<<
        gatk CollectAlignmentSummaryMetrics \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --INPUT ~{bam} \
        --OUTPUT ~{prefix}.alignment_summary_metrics.txt \
        --REFERENCE_SEQUENCE ~{reference} \
        --VALIDATION_STRINGENCY LENIENT
    >>>

    output {
        File alignment_summary_metrics = "~{prefix}.alignment_summary_metrics.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

task CollectInsertSizeMetrics {
    input {
        File bam
        File bai

        Int cpu = 2
        Int memory_gb = 16
        Int disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
    }

    String prefix = basename(bam, ".sorted.bam")

    command <<<
        gatk CollectInsertSizeMetrics \
        --INPUT ~{bam} \
        --OUTPUT ~{prefix}.insert_size_metrics.txt \
        --Histogram_FILE ~{prefix}.insert_size_plot.pdf \
        --VALIDATION_STRINGENCY LENIENT
    >>>

    output {
        File insert_size_metrics = "~{prefix}.insert_size_metrics.txt"
        File insert_size_plot = "~{prefix}.insert_size_plot.pdf"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

task CollectUMIDuplicationMetrics {
    input {
        String output_basename
        File umi_group_data

        Int cpu = 1
        Int memory_gb = 8
        Int disk_size_gb = 32
    }

    command <<<
        set -e
        python3 <<CODE

        umi_group_data_dict = {}

        # Tab-separated list of columns: [family_size, count, fraction, fraction_gt_or_eq_family_size]
        with open("~{umi_group_data}", 'r') as infile:
            header = infile.readline()
            for line in infile:
                line = line.rstrip()
                columns = line.split('\t')
                umi_group_data_dict[int(columns[0])] = (int(columns[1]), float(columns[2]), float(columns[3]))

        num_fragments_total = 0
        num_fragments_unique = 0
        for family_size, value_tuple in umi_group_data_dict.items():
            num_fragments_total = num_fragments_total + (family_size * value_tuple[0])
            num_fragments_unique = num_fragments_unique + value_tuple[0]

        percent_duplication = 100 * (1 - (num_fragments_unique / num_fragments_total))

        with open("~{output_basename}.umi_duplication_metrics.tsv", 'w') as f:
            f.write("PERCENT_DUPLICATION" + "\t" + str(percent_duplication) + "\n")
            f.write("ESTIMATED_LIBRARY_SIZE" + "\t" + str(num_fragments_unique))
        CODE
    >>>

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    }

    output {
        File umi_duplication_metrics = "~{output_basename}.umi_duplication_metrics.tsv"
    }
}

task CollectHsMetrics {
    input {
        File bam
        File bai
        File reference
        File reference_fai
        File reference_dict
        File bait_interval_list
        File target_interval_list
        String bait_set_name
        String output_prefix

        Int cpu = 2
        Int memory_gb = 32
        Int disk_size_gb = 512
    }

    command <<<
        gatk CollectHsMetrics \
        --BAIT_SET_NAME ~{bait_set_name} \
        --BAIT_INTERVALS ~{bait_interval_list} \
        --TARGET_INTERVALS ~{target_interval_list} \
        --INPUT ~{bam} \
        --OUTPUT ~{output_prefix}.hs_metrics.txt \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --REFERENCE_SEQUENCE ~{reference} \
        --COVERAGE_CAP 100000 \
        --PER_TARGET_COVERAGE ~{output_prefix}.per_target_coverage.txt \
        --VALIDATION_STRINGENCY LENIENT
    >>>

    output {
        File hs_metrics = "~{output_prefix}.hs_metrics.txt"
        File per_target_coverage = "~{output_prefix}.per_target_coverage.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

task CollectDuplexSeqMetrics {
    input {
        File bam

        Int cpu = 2
        Int memory_gb = 16
        Int disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
    }

    String prefix = basename(bam, ".duplex.umi_grouped.bam")

    command <<<
        fgbio CollectDuplexSeqMetrics \
        --input ~{bam} \
        --output ~{prefix}
    >>>

    output {
        File family_sizes = "~{prefix}.family_sizes.txt"
        File duplex_family_sizes = "~{prefix}.duplex_family_sizes.txt"
        File duplex_yield_metrics = "~{prefix}.duplex_yield_metrics.txt"
        File umi_counts = "~{prefix}.umi_counts.txt"
        File duplex_qc = "~{prefix}.duplex_qc.pdf"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/fgbio@sha256:b6869a0ae243d9f1b183e4a986fbe0853df2a56a1c6d7c0fec2965b6d8a7af1d"
    }
}

workflow HPVDeepSeekGenotyping {
    input {
        String output_basename
        File r1_fastq
        File r2_fastq
        File human_snp_targets_bed
        File reference
        File reference_fai
        File reference_dict
        File bwa_idx_amb
        File bwa_idx_ann
        File bwa_idx_bwt
        File bwa_idx_pac
        File bwa_idx_sa
        File hpv_bait_interval_list
        File hpv_target_interval_list
        File hg38_bait_interval_list
        File hg38_target_interval_list
        File low_risk_hpv_genotypes
        String bait_set_name
        String read_group_id
        String read_group_sample_name
        String read_group_library_name = "LB_TEST"
        String read_group_platform = "ILLUMINA"
        String read_group_platform_unit = "PU_TEST"
        String read_group_description = "KAPA_TE"
        String read_structure
    }

    call FastQC as PreTrimmedFastQC {
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq
    }

    call FastqToUbam {
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            output_basename = output_basename,
            read_group_id = read_group_id,
            read_group_sample_name = read_group_sample_name,
            read_group_library_name = read_group_library_name,
            read_group_platform = read_group_platform,
            read_group_platform_unit = read_group_platform_unit,
            read_group_description = read_group_description
    }

    call ExtractUMIs {
        input:
            input_ubam = FastqToUbam.ubam,
            read_structure = read_structure,
            output_basename = output_basename
    }

    call UmiExtractedBamToFastq {
        input:
            umi_extracted_bam = ExtractUMIs.umi_extracted_bam,
            output_basename = output_basename
    }

    call TrimAndFilter {
        input:
            fastq1 = UmiExtractedBamToFastq.umi_extracted_fastq1,
            fastq2 = UmiExtractedBamToFastq.umi_extracted_fastq2,
            output_basename = output_basename
    }

    call FastQC as PostTrimmedFastQC {
        input:
            r1_fastq = TrimAndFilter.fastq1_trimmed,
            r2_fastq = TrimAndFilter.fastq2_trimmed
    }

    call BwaMem as AlignReads {
        input:
            fastq1 = TrimAndFilter.fastq1_trimmed,
            fastq2 = TrimAndFilter.fastq2_trimmed,
            reference = reference,
            bwa_idx_amb = bwa_idx_amb,
            bwa_idx_ann = bwa_idx_ann,
            bwa_idx_bwt = bwa_idx_bwt,
            bwa_idx_pac = bwa_idx_pac,
            bwa_idx_sa = bwa_idx_sa,
            read_group_id = read_group_id,
            read_group_sample_name = read_group_sample_name,
            read_group_library_name = read_group_library_name,
            read_group_platform = read_group_platform,
            read_group_platform_unit = read_group_platform_unit,
            read_group_description = read_group_description,
            output_basename = output_basename
    }

    call SortAndIndexBam {
        input:
            bam = AlignReads.bam
    }

    call CollectAlignmentSummaryMetrics as PreConsensusAlignmentSummaryMetrics {
        input:
            bam = SortAndIndexBam.sorted_bam,
            bai = SortAndIndexBam.sorted_bam_index,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict
    }

    call CollectInsertSizeMetrics as PreConsensusInsertSizeMetrics {
        input:
            bam = SortAndIndexBam.sorted_bam,
            bai = SortAndIndexBam.sorted_bam_index
    }

    call CollectHsMetrics as CollectHsMetricsRawHPV {
        input:
            bam = SortAndIndexBam.sorted_bam,
            bai = SortAndIndexBam.sorted_bam_index,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bait_interval_list = hpv_bait_interval_list,
            target_interval_list = hpv_target_interval_list,
            bait_set_name = bait_set_name,
            output_prefix = output_basename + ".raw.hpv"
    }

    call CollectHsMetrics as CollectHsMetricsRawHg38 {
        input:
            bam = SortAndIndexBam.sorted_bam,
            bai = SortAndIndexBam.sorted_bam_index,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bait_interval_list = hg38_bait_interval_list,
            target_interval_list = hg38_target_interval_list,
            bait_set_name = bait_set_name,
            output_prefix = output_basename + ".raw.hg38"
    }

    call MergeBamAlignment as MergeBAMs {
        input:
            aligned_bam = SortAndIndexBam.sorted_bam,
            unmapped_bam = ExtractUMIs.umi_extracted_bam,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            output_basename = output_basename + ".merged",
            extra_args = "--SORT_ORDER queryname --ALIGNED_READS_ONLY true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ALIGNER_PROPER_PAIR_FLAGS true --CLIP_OVERLAPPING_READS false --ATTRIBUTES_TO_RETAIN X0 --ATTRIBUTES_TO_REMOVE NM --ATTRIBUTES_TO_REMOVE MD"
    }

    call FilterAndGroupReadsByUMI as FilterAndGroupReadsByUMISimplex {
        input:
            merged_bam = MergeBAMs.merged_bam,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            is_duplex = false,
            output_basename = output_basename
    }

    call FilterAndGroupReadsByUMI as FilterAndGroupReadsByUMIDuplex {
        input:
            merged_bam = MergeBAMs.merged_bam,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            is_duplex = true,
            output_basename = output_basename
    }

    call CollectUMIDuplicationMetrics as CollectUMIDuplicationMetricsSimplex {
        input:
            umi_group_data = FilterAndGroupReadsByUMISimplex.umi_group_data,
            output_basename = output_basename + ".simplex"
    }

    call CollectUMIDuplicationMetrics as CollectUMIDuplicationMetricsDuplex {
        input:
            umi_group_data = FilterAndGroupReadsByUMIDuplex.umi_group_data,
            output_basename = output_basename + ".duplex"
    }

    call CallMolecularConsensusReads {
        input:
            umi_grouped_bam = FilterAndGroupReadsByUMISimplex.umi_grouped_bam,
            read_group_id = read_group_id,
            output_basename = output_basename
    }

    call CallDuplexConsensusReads {
        input:
            umi_grouped_bam = FilterAndGroupReadsByUMIDuplex.umi_grouped_bam,
            read_group_id = read_group_id,
            output_basename = output_basename
    }

    call CollectDuplexSeqMetrics {
        input:
            bam = FilterAndGroupReadsByUMIDuplex.umi_grouped_bam
    }

    call ConsensusBamToFastq as SimplexConsensusBamToFastq {
        input:
            umi_consensus_unmapped_bam = CallMolecularConsensusReads.umi_consensus_unmapped_bam,
            output_basename = output_basename
    }

    call ConsensusBamToFastq as DuplexConsensusBamToFastq {
        input:
            umi_consensus_unmapped_bam = CallDuplexConsensusReads.umi_consensus_unmapped_bam,
            output_basename = output_basename
    }

    # Align consensus reads to the reference genome:
    # bwa mem -t THREADS \
    # -R '@RG\tID:A\tDS:KAPA_TE\tPL:ILLUMINA\tLB:lib1\tSM:sample1\tPU:unit1' \
    # -v 3 -Y -M -K 100000000 \
    # reference.fasta \
    # consensus_unmapped_R1.fastq consensus_unmapped_R2.fastq | \
    # samtools view -bh - > consensus_mapped_unsorted.bam
    call BwaMem as AlignSimplexConsensusReads {
        input:
            fastq1 = SimplexConsensusBamToFastq.consensus_unmapped_fastq1,
            fastq2 = SimplexConsensusBamToFastq.consensus_unmapped_fastq2,
            reference = reference,
            bwa_idx_amb = bwa_idx_amb,
            bwa_idx_ann = bwa_idx_ann,
            bwa_idx_bwt = bwa_idx_bwt,
            bwa_idx_pac = bwa_idx_pac,
            bwa_idx_sa = bwa_idx_sa,
            read_group_id = read_group_id,
            read_group_sample_name = read_group_sample_name,
            read_group_library_name = read_group_library_name,
            read_group_platform = read_group_platform,
            read_group_platform_unit = read_group_platform_unit,
            read_group_description = read_group_description,
            soft_clip_supplementary_alignments = true,
            output_basename = output_basename + ".simplex.consensus"
    }

    call BwaMem as AlignDuplexConsensusReads {
        input:
            fastq1 = DuplexConsensusBamToFastq.consensus_unmapped_fastq1,
            fastq2 = DuplexConsensusBamToFastq.consensus_unmapped_fastq2,
            reference = reference,
            bwa_idx_amb = bwa_idx_amb,
            bwa_idx_ann = bwa_idx_ann,
            bwa_idx_bwt = bwa_idx_bwt,
            bwa_idx_pac = bwa_idx_pac,
            bwa_idx_sa = bwa_idx_sa,
            read_group_id = read_group_id,
            read_group_sample_name = read_group_sample_name,
            read_group_library_name = read_group_library_name,
            read_group_platform = read_group_platform,
            read_group_platform_unit = read_group_platform_unit,
            read_group_description = read_group_description,
            soft_clip_supplementary_alignments = true,
            output_basename = output_basename + ".duplex.consensus"
    }

    call GATKSortBam as GATKSortBamSimplexConsensusAligned {
        input:
            bam = AlignSimplexConsensusReads.bam
    }

    call GATKSortBam as GATKSortBamDuplexConsensusAligned {
        input:
            bam = AlignDuplexConsensusReads.bam
    }

    call GATKSortBam as GATKSortBamSimplexConsensusUnmapped {
        input:
            bam = CallMolecularConsensusReads.umi_consensus_unmapped_bam
    }

    call GATKSortBam as GATKSortBamDuplexConsensusUnmapped {
        input:
            bam = CallDuplexConsensusReads.umi_consensus_unmapped_bam
    }

    call MergeBamAlignment as MergeConsensusSimplex {
        input:
            aligned_bam = GATKSortBamSimplexConsensusAligned.sorted_bam,
            unmapped_bam = GATKSortBamSimplexConsensusUnmapped.sorted_bam,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            output_basename = output_basename + ".simplex",
            extra_args = "--SORT_ORDER coordinate --ATTRIBUTES_TO_RETAIN X0 --ATTRIBUTES_TO_RETAIN RX --ADD_MATE_CIGAR true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ALIGNER_PROPER_PAIR_FLAGS true --CLIP_OVERLAPPING_READS false"
    }

    call MergeBamAlignment as MergeConsensusDuplex {
        input:
            aligned_bam = GATKSortBamDuplexConsensusAligned.sorted_bam,
            unmapped_bam = GATKSortBamDuplexConsensusUnmapped.sorted_bam,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            output_basename = output_basename + ".duplex",
            extra_args = "--SORT_ORDER coordinate --ATTRIBUTES_TO_RETAIN X0 --ATTRIBUTES_TO_RETAIN RX --ADD_MATE_CIGAR true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ALIGNER_PROPER_PAIR_FLAGS true --CLIP_OVERLAPPING_READS false"
    }

    call SortAndIndexBam as SortAndIndexSimplexBam {
        input:
            bam = MergeConsensusSimplex.merged_bam
    }

    call SortAndIndexBam as SortAndIndexDuplexBam {
        input:
            bam = MergeConsensusDuplex.merged_bam
    }

    call CollectAlignmentSummaryMetrics as PostConsensusAlignmentSummaryMetrics {
        input:
            bam = SortAndIndexSimplexBam.sorted_bam,
            bai = SortAndIndexSimplexBam.sorted_bam_index,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict
    }

    call CollectInsertSizeMetrics as PostConsensusInsertSizeMetrics {
        input:
            bam = SortAndIndexSimplexBam.sorted_bam,
            bai = SortAndIndexSimplexBam.sorted_bam_index
    }

    call CollectHsMetrics as CollectHsMetricsSimplexHPV {
        input:
            bam = SortAndIndexSimplexBam.sorted_bam,
            bai = SortAndIndexSimplexBam.sorted_bam_index,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bait_interval_list = hpv_bait_interval_list,
            target_interval_list = hpv_target_interval_list,
            bait_set_name = bait_set_name,
            output_prefix = output_basename + ".simplex.hpv"
    }

    call CollectHsMetrics as CollectHsMetricsDuplexHPV {
        input:
            bam = SortAndIndexDuplexBam.sorted_bam,
            bai = SortAndIndexDuplexBam.sorted_bam_index,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bait_interval_list = hpv_bait_interval_list,
            target_interval_list = hpv_target_interval_list,
            bait_set_name = bait_set_name,
            output_prefix = output_basename + ".duplex.hpv"
    }

    call CollectHsMetrics as CollectHsMetricsSimplexHg38 {
        input:
            bam = SortAndIndexSimplexBam.sorted_bam,
            bai = SortAndIndexSimplexBam.sorted_bam_index,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bait_interval_list = hg38_bait_interval_list,
            target_interval_list = hg38_target_interval_list,
            bait_set_name = bait_set_name,
            output_prefix = output_basename + ".simplex.hg38"
    }

    call CollectHsMetrics as CollectHsMetricsDuplexHg38 {
        input:
            bam = SortAndIndexDuplexBam.sorted_bam,
            bai = SortAndIndexDuplexBam.sorted_bam_index,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bait_interval_list = hg38_bait_interval_list,
            target_interval_list = hg38_target_interval_list,
            bait_set_name = bait_set_name,
            output_prefix = output_basename + ".duplex.hg38"
    }

    call SamtoolsCoverage {
        input:
            bam = SortAndIndexDuplexBam.sorted_bam,
            output_basename = output_basename
    }

    call DetermineHPVStatus {
        input:
            coverage = SamtoolsCoverage.coverage,
            low_risk_hpv_genotypes = low_risk_hpv_genotypes,
            output_basename = output_basename
    }

    call GenotypeSNPsHuman {
        input:
            bam = SortAndIndexSimplexBam.sorted_bam,
            bai = SortAndIndexSimplexBam.sorted_bam_index,
            human_snp_targets_bed = human_snp_targets_bed,
            reference = reference,
            output_basename = output_basename
    }

    output {
        File raw_bam = SortAndIndexBam.sorted_bam
        File raw_bam_index = SortAndIndexBam.sorted_bam_index
        File simplex_bam = SortAndIndexSimplexBam.sorted_bam
        File simplex_bam_index = SortAndIndexSimplexBam.sorted_bam_index
        File duplex_bam = SortAndIndexDuplexBam.sorted_bam
        File duplex_bam_index = SortAndIndexDuplexBam.sorted_bam_index
        File simplex_umi_grouped_bam = FilterAndGroupReadsByUMISimplex.umi_grouped_bam
        File simplex_umi_group_data = FilterAndGroupReadsByUMISimplex.umi_group_data
        File duplex_umi_grouped_bam = FilterAndGroupReadsByUMIDuplex.umi_grouped_bam
        File duplex_umi_group_data = FilterAndGroupReadsByUMIDuplex.umi_group_data
        File simplex_umi_duplication_metrics = CollectUMIDuplicationMetricsSimplex.umi_duplication_metrics
        File duplex_umi_duplication_metrics = CollectUMIDuplicationMetricsDuplex.umi_duplication_metrics
        File vcf = GenotypeSNPsHuman.vcf
        File coverage = SamtoolsCoverage.coverage
        File hpv_status = DetermineHPVStatus.hpv_status
        File fastp_report_html = TrimAndFilter.fastp_report_html
        File fastp_report_json = TrimAndFilter.fastp_report_json
        File pre_trimmed_r1_fastqc_html = PreTrimmedFastQC.r1_fastqc_html
        File pre_trimmed_r2_fastqc_html = PreTrimmedFastQC.r2_fastqc_html
        File post_trimmed_r1_fastqc_html = PostTrimmedFastQC.r1_fastqc_html
        File post_trimmed_r2_fastqc_html = PostTrimmedFastQC.r2_fastqc_html
        File pre_consensus_alignment_summary_metrics = PreConsensusAlignmentSummaryMetrics.alignment_summary_metrics
        File pre_consensus_insert_size_metrics = PreConsensusInsertSizeMetrics.insert_size_metrics
        File pre_consensus_insert_size_plot = PreConsensusInsertSizeMetrics.insert_size_plot
        File post_consensus_alignment_summary_metrics = PostConsensusAlignmentSummaryMetrics.alignment_summary_metrics
        File post_consensus_insert_size_metrics = PostConsensusInsertSizeMetrics.insert_size_metrics
        File post_consensus_insert_size_plot = PostConsensusInsertSizeMetrics.insert_size_plot
        File raw_hpv_hs_metrics = CollectHsMetricsRawHPV.hs_metrics
        File raw_hpv_per_target_coverage = CollectHsMetricsRawHPV.per_target_coverage
        File raw_hg38_hs_metrics = CollectHsMetricsRawHg38.hs_metrics
        File raw_hg38_per_target_coverage = CollectHsMetricsRawHg38.per_target_coverage
        File simplex_hpv_hs_metrics = CollectHsMetricsSimplexHPV.hs_metrics
        File simplex_hpv_per_target_coverage = CollectHsMetricsSimplexHPV.per_target_coverage
        File simplex_hg38_hs_metrics = CollectHsMetricsSimplexHg38.hs_metrics
        File simplex_hg38_per_target_coverage = CollectHsMetricsSimplexHg38.per_target_coverage
        File duplex_hpv_hs_metrics = CollectHsMetricsDuplexHPV.hs_metrics
        File duplex_hpv_per_target_coverage = CollectHsMetricsDuplexHPV.per_target_coverage
        File duplex_hg38_hs_metrics = CollectHsMetricsDuplexHg38.hs_metrics
        File duplex_hg38_per_target_coverage = CollectHsMetricsDuplexHg38.per_target_coverage
        File family_sizes = CollectDuplexSeqMetrics.family_sizes
        File duplex_family_sizes = CollectDuplexSeqMetrics.duplex_family_sizes
        File duplex_yield_metrics = CollectDuplexSeqMetrics.duplex_yield_metrics
        File umi_counts = CollectDuplexSeqMetrics.umi_counts
        File duplex_qc = CollectDuplexSeqMetrics.duplex_qc
    }
}