version 1.0

task FastQC {
    input {
        File r1_fastq
        File r2_fastq
        Int? fastqc_thread_memory = 4096
        Int? cpu = 2
        Int? num_threads = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2 * (size(r1_fastq, "GiB") + size(r2_fastq, "GiB"))) + 50)
        Int? min_ssd_size_gb = 512
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
        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * (size(r1_fastq, "GiB") + size(r2_fastq, "GiB"))) + 50)
        Int? min_ssd_size_gb = 512
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

# Extract the UMI sequence from the first 3 bases of each read, skip the next 3 bases, append a 'T',
# and add the resulting UMI to the RX tag and read name in the output BAM.
#
# --molecular-index-tags RX: Tells fgbio to place the extracted UMI into the RX tag of each read.
# --read-structure 3M3S+T 3M3S+T: Specifies the regular expressions used to extract the UMI from each read sequence:
# 3M3S+T means:
# 3M: Match the first 3 bases (UMI).
# 3S: Skip (soft clip) the next 3 bases.
# +T: Append this to the UMI (the sequence of "T").
# This regex is applied to both R1 and R2 (paired-end reads).
# --annotate-read-names true: Indicates that the UMI should be appended to the read name (QNAME) in addition to the RX tag.
task ExtractUMIs {
    input {
        String output_basename
        File input_ubam
        String read_group_tag = "RX"
        String read_structure = "3M3S+T"
        String append_umi_to_qname = "true"
        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(input_ubam, "GiB")) + 50)
        Int? min_ssd_size_gb = 512
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
        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(umi_extracted_bam, "GiB")) + 50)
        Int? min_ssd_size_gb = 512
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
        Int? cpu = 3
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * (size(fastq1, "GiB") + size(fastq2, "GiB"))) + 50)
        Int? min_ssd_size_gb = 512
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
        File bwa_idx_amb
        File bwa_idx_ann
        File bwa_idx_bwt
        File bwa_idx_pac
        File bwa_idx_sa
        String read_group_id
        String read_group_sample_name
        String read_group_library_name
        String read_group_platform
        String read_group_platform_unit
        String read_group_description
        Boolean soft_clip_supplementary_alignments = false
        Int? cpu = 32
        Int? num_threads = 32
        Int? memory_gb = 64
        Int? disk_size_gb = ceil((3 * (size(fastq1, "GiB") + size(fastq2, "GiB"))) + size(reference, "GiB") + 100)
        Int? min_ssd_size_gb = 512
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
        String? samtools_thread_memory = "1024M"
        Int? cpu = 7
        Int? num_threads = 7
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
        Int? min_ssd_size_gb = 512
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
        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
        Int? min_ssd_size_gb = 512
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

# 1. Merge BWA-aligned reads with original UMI-tagged reads to produce a BAM file that has alignment info + UMI tags, necessary for UMI-aware deduplication.
# 2. Filter reads with MAPQ < 1 and reads not mapped in proper pair
# 3. Group together reads that share the same UMI sequence (or a similar one within a defined edit distance) and originate from the same genomic location—used for UMI-aware deduplication.
# fgbio:
#--input: Input BAM file containing reads aligned to the reference with UMI tags (RX) and filtered.
#--output: Output BAM where each read is tagged with a group ID (MI tag) indicating its molecular family.
#--strategy=adjacency: Uses an adjacency graph to group UMIs—UMIs within one base mismatch (--edits=1) are grouped if one is more abundant.
#--edits=1: Allows grouping of UMIs within 1 base mismatch (to account for sequencing errors).
#-t RX: Specifies the tag where UMIs are stored (here, the standard RX tag).
#-f ...umi_group_data.txt: Outputs grouping statistics and metadata for each molecular family.
task MergeBAMsAndGroupUMIs {
    input {
        String output_basename
        File aligned_bam
        File unmapped_umi_extracted_bam
        File reference
        File reference_fai
        File reference_dict
        Boolean call_duplex_consensus

        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(aligned_bam, "GiB") + size(unmapped_umi_extracted_bam, "GiB")) + 100)
        Int? min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    String strategy = if call_duplex_consensus then "paired" else "adjacency"

    command <<<
        gatk MergeBamAlignment \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ATTRIBUTES_TO_REMOVE NM \
        --ATTRIBUTES_TO_REMOVE MD \
        --ALIGNED_BAM ~{aligned_bam} \
        --UNMAPPED_BAM ~{unmapped_umi_extracted_bam} \
        --OUTPUT ~{output_basename}.merged.bam \
        --REFERENCE_SEQUENCE ~{reference} \
        --SORT_ORDER queryname \
        --ALIGNED_READS_ONLY true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --CLIP_OVERLAPPING_READS false

        samtools view -f 2 -q 1 -bh ~{output_basename}.merged.bam > ~{output_basename}.merged.filtered.bam

        fgbio GroupReadsByUmi \
        --input ~{output_basename}.merged.filtered.bam \
        --output ~{output_basename}.umi_grouped.bam \
        --strategy ~{strategy} \
        --edits 1 \
        --raw-tag RX \
        --family-size-histogram ~{output_basename}.umi_group_data.txt
    >>>

    output {
        File umi_grouped_bam = "~{output_basename}.umi_grouped.bam"
        File umi_group_data = "~{output_basename}.umi_group_data.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# Merge aligned reads with UMI and original tags:
# gatk MergeBamAlignment \
# --ALIGNED_BAM consensus_mapped.bam \
# --UNMAPPED_BAM consensus_unmapped_sorted.bam \
# --OUTPUT deduped.bam \
# --REFERENCE_SEQUENCE reference.fasta \
# --SORT_ORDER coordinate \
# --ATTRIBUTES_TO_RETAIN X0 \
# --ATTRIBUTES_TO_RETAIN RX \
# --ADD_MATE_CIGAR true \
# --MAX_INSERTIONS_OR_DELETIONS -1 \
# --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
# --ALIGNER_PROPER_PAIR_FLAGS true \
# --CLIP_OVERLAPPING_READS false
task MergeConsensus {
    input {
        String output_basename
        File consensus_aligned_bam
        File consensus_unmapped_bam
        File reference
        File reference_fai
        File reference_dict
        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(consensus_aligned_bam, "GiB") + size(consensus_unmapped_bam, "GiB")) + 100)
        Int? min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    command <<<
        gatk MergeBamAlignment \
        --ALIGNED_BAM ~{consensus_aligned_bam} \
        --UNMAPPED_BAM ~{consensus_unmapped_bam} \
        --OUTPUT ~{output_basename}.deduped.bam \
        --REFERENCE_SEQUENCE ~{reference} \
        --SORT_ORDER coordinate \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ATTRIBUTES_TO_RETAIN RX \
        --ADD_MATE_CIGAR true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --CLIP_OVERLAPPING_READS false
    >>>

    output {
        File deduped_bam = "~{output_basename}.deduped.bam"
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
        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(umi_grouped_bam, "GiB")) + 100)
    }

    command <<<
        fgbio CallMolecularConsensusReads \
        --input ~{umi_grouped_bam} \
        --output ~{output_basename}.umi_consensus.unmapped.bam \
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
        File umi_consensus_unmapped_bam = "~{output_basename}.umi_consensus.unmapped.bam"
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
        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(umi_grouped_bam, "GiB")) + 100)
    }

    command <<<
        fgbio CallDuplexConsensusReads \
        --input ~{umi_grouped_bam} \
        --output ~{output_basename}.umi_consensus.unmapped.bam \
        --error-rate-post-umi 40 \
        --error-rate-pre-umi 45 \
        --min-reads 1 \
        --max-reads-per-strand 50 \
        --min-input-base-quality 20 \
        --read-name-prefix 'consensus' \
        --read-group-id ~{read_group_id} \
        --consensus-call-overlapping-bases true \
        --stats ~{output_basename}.consensus.stats.txt
    >>>

    output {
        File umi_consensus_unmapped_bam = "~{output_basename}.umi_consensus.unmapped.bam"
        File consensus_stats = "~{output_basename}.consensus.stats.txt"
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
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(umi_consensus_unmapped_bam, "GiB")) + 50)
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
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2 * size(bam, "GiB")) + 50)
    }

    command <<<
        samtools coverage ~{bam} | awk '{if(NR==1){printf "%s\t%s\t%s\n",$1,$4,$6} else {printf "%s\t%s\t%s\n",$1,$4,$6}}' > ~{output_basename}.coverage.txt
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
# Sample is considered HPV+ if both thresholds are met:
# Read count ≥ 10
# Percentage of HPV genome covered by read alignment ≥ 10%
task DetermineHPVStatus {
    input {
        File coverage
        Int? cpu = 1
        Int? memory_gb = 8
        Int? disk_size_gb = 32
    }

    command <<<
        set -e
        python3 <<CODE

        coverage_dict = {}
        with open("~{coverage}", 'r') as infile:
            header = infile.readline()
            for line in infile:
                line = line.rstrip()
                columns = line.split('\t')
                if columns[0].startswith("HPV"):
                    coverage_dict[columns[0]] = (int(columns[1]), float(columns[2]))

        coverage_sorted = sorted(coverage_dict.items(), key = lambda item: item[1], reverse = True)
        max_elem = coverage_sorted[0]

        with open("top_hpv_contig.txt", 'w') as f:
            f.write(max_elem[0])

        with open("top_hpv_num_reads.txt", 'w') as f:
            f.write(str(max_elem[1][0]))

        with open("top_hpv_coverage.txt", 'w') as f:
            f.write(str(max_elem[1][1]))

        with open("is_hpv_positive.txt", 'w') as f:
            if max_elem[1][0] >= 10 and max_elem[1][1] >= 10.0:
                f.write("true")
            else:
                f.write("false")

        with open("secondary_hpv_types.txt", 'w') as f:
            output_string = ""
            for i in range(1, len(coverage_sorted)):
                elem = coverage_sorted[i]
                if elem[1][0] >= 10 and elem[1][1] >= 10.0:
                    output_string = output_string + str(elem[0]) + ":" + str(elem[1][0]) + ":" + str(elem[1][1]) + ","

            if len(output_string) > 0:
                output_string = output_string[:-1]
            f.write(output_string)
        CODE
    >>>

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    }

    output {
        String top_hpv_contig = read_string("top_hpv_contig.txt")
        Int top_hpv_num_reads = read_int("top_hpv_num_reads.txt")
        Float top_hpv_coverage = read_float("top_hpv_coverage.txt")
        Boolean is_hpv_positive = read_boolean("is_hpv_positive.txt")
        String secondary_hpv_types = read_string("secondary_hpv_types.txt")
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
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(bam, "GiB")) + 50)
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
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
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

task Flagstat {
    input {
        File bam
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
    }

    String prefix = basename(bam, ".sorted.bam")

    command <<<
        samtools flagstat ~{bam} > ~{prefix}.flagstat.txt
    >>>

    output {
        File flagstat = "~{prefix}.flagstat.txt"
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
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
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

task CountOnTargetReads {
    input {
        File bam
        File bai
        File reference
        File reference_fai
        File reference_dict
        File capture_targets_bed
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
    }

    String prefix = basename(bam, ".sorted.bam")

    command <<<
        gatk CountReads \
        --reference ~{reference} \
        --input ~{bam} \
        --intervals ~{capture_targets_bed} \
        --read-filter MappedReadFilter \
        --read-filter NotSecondaryAlignmentReadFilter > ~{prefix}.ontarget_reads.txt
    >>>

    output {
        File ontarget_reads = "~{prefix}.ontarget_reads.txt"
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
        Int? cpu = 1
        Int? memory_gb = 8
        Int? disk_size_gb = 32
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

task CollectHybridSelectionMetrics {
    input {
        File bam
        File bai
        File reference
        File reference_fai
        File reference_dict
        File bait_interval_list
        File target_interval_list
        String bait_set_name
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
    }

    String prefix = basename(bam, ".sorted.bam")

    command <<<
        gatk CollectHsMetrics \
        --BAIT_INTERVALS ~{bait_interval_list} \
        --BAIT_SET_NAME ~{bait_set_name} \
        --TARGET_INTERVALS ~{target_interval_list} \
        --INPUT ~{bam} \
        --OUTPUT ~{prefix}.hs_metrics.txt \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --REFERENCE_SEQUENCE ~{reference} \
        --COVERAGE_CAP 100000 \
        --PER_BASE_COVERAGE ~{prefix}.per_base_coverage.txt \
        --VALIDATION_STRINGENCY LENIENT
    >>>

    output {
        File hs_metrics = "~{prefix}.hs_metrics.txt"
        File per_base_coverage = "~{prefix}.per_base_coverage.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
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
        File capture_targets_bed
        File bait_interval_list
        File target_interval_list
        String bait_set_name
        String read_group_id
        String read_group_sample_name
        String read_group_library_name = "LB_TEST"
        String read_group_platform = "ILLUMINA"
        String read_group_platform_unit = "PU_TEST"
        String read_group_description = "KAPA_TE"
        Boolean call_duplex_consensus = false
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

    call Flagstat as PreConsensusFlagstat {
        input:
            bam = SortAndIndexBam.sorted_bam
    }

    call CollectInsertSizeMetrics as PreConsensusInsertSizeMetrics {
        input:
            bam = SortAndIndexBam.sorted_bam,
            bai = SortAndIndexBam.sorted_bam_index
    }

    call CountOnTargetReads as PreConsensusCountOnTargetReads {
        input:
            bam = SortAndIndexBam.sorted_bam,
            bai = SortAndIndexBam.sorted_bam_index,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            capture_targets_bed = capture_targets_bed
    }

    call CollectHybridSelectionMetrics as PreConsensusHybridSelectionMetrics {
        input:
            bam = SortAndIndexBam.sorted_bam,
            bai = SortAndIndexBam.sorted_bam_index,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bait_interval_list = bait_interval_list,
            target_interval_list = target_interval_list,
            bait_set_name = bait_set_name
    }

    call MergeBAMsAndGroupUMIs {
        input:
            aligned_bam = SortAndIndexBam.sorted_bam,
            unmapped_umi_extracted_bam = ExtractUMIs.umi_extracted_bam,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            call_duplex_consensus = call_duplex_consensus,
            output_basename = output_basename
    }

    call CollectUMIDuplicationMetrics {
        input:
            umi_group_data = MergeBAMsAndGroupUMIs.umi_group_data,
            output_basename = output_basename
    }

    if(!call_duplex_consensus) {
        call CallMolecularConsensusReads {
            input:
                umi_grouped_bam = MergeBAMsAndGroupUMIs.umi_grouped_bam,
                read_group_id = read_group_id,
                output_basename = output_basename
        }
    }

    if(call_duplex_consensus) {
        call CallDuplexConsensusReads {
            input:
                umi_grouped_bam = MergeBAMsAndGroupUMIs.umi_grouped_bam,
                read_group_id = read_group_id,
                output_basename = output_basename
        }
    }

    File umi_consensus_unmapped_bam = select_first([CallDuplexConsensusReads.umi_consensus_unmapped_bam, CallMolecularConsensusReads.umi_consensus_unmapped_bam])

    call ConsensusBamToFastq {
        input:
            umi_consensus_unmapped_bam = umi_consensus_unmapped_bam,
            output_basename = output_basename
    }

    # Align consensus reads to the reference genome:
    # bwa mem -t THREADS \
    # -R '@RG\tID:A\tDS:KAPA_TE\tPL:ILLUMINA\tLB:lib1\tSM:sample1\tPU:unit1' \
    # -v 3 -Y -M -K 100000000 \
    # reference.fasta \
    # consensus_unmapped_R1.fastq consensus_unmapped_R2.fastq | \
    # samtools view -bh - > consensus_mapped_unsorted.bam
    String consensus_basename = output_basename + ".consensus"
    call BwaMem as AlignConsensusReads {
        input:
            fastq1 = ConsensusBamToFastq.consensus_unmapped_fastq1,
            fastq2 = ConsensusBamToFastq.consensus_unmapped_fastq2,
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
            output_basename = consensus_basename
    }

    call GATKSortBam as GATKSortBamConsensusAligned {
        input:
            bam = AlignConsensusReads.bam
    }

    call GATKSortBam as GATKSortBamConsensusUnmapped{
        input:
            bam = umi_consensus_unmapped_bam
    }

    call MergeConsensus {
        input:
            consensus_aligned_bam = GATKSortBamConsensusAligned.sorted_bam,
            consensus_unmapped_bam = GATKSortBamConsensusUnmapped.sorted_bam,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            output_basename = output_basename
    }

    call SortAndIndexBam as SortAndIndexFinalBam {
        input:
            bam = MergeConsensus.deduped_bam
    }

    call CollectAlignmentSummaryMetrics as PostConsensusAlignmentSummaryMetrics {
        input:
            bam = SortAndIndexFinalBam.sorted_bam,
            bai = SortAndIndexFinalBam.sorted_bam_index,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict
    }

    call Flagstat as PostConsensusFlagstat {
        input:
            bam = SortAndIndexFinalBam.sorted_bam
    }

    call CollectInsertSizeMetrics as PostConsensusInsertSizeMetrics {
        input:
            bam = SortAndIndexFinalBam.sorted_bam,
            bai = SortAndIndexFinalBam.sorted_bam_index
    }

    call CountOnTargetReads as PostConsensusCountOnTargetReads {
        input:
            bam = SortAndIndexFinalBam.sorted_bam,
            bai = SortAndIndexFinalBam.sorted_bam_index,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            capture_targets_bed = capture_targets_bed
    }

    call CollectHybridSelectionMetrics as PostConsensusHybridSelectionMetrics {
        input:
            bam = SortAndIndexFinalBam.sorted_bam,
            bai = SortAndIndexFinalBam.sorted_bam_index,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bait_interval_list = bait_interval_list,
            target_interval_list = target_interval_list,
            bait_set_name = bait_set_name
    }

    call SamtoolsCoverage {
        input:
            bam = SortAndIndexFinalBam.sorted_bam,
            output_basename = output_basename
    }

    call DetermineHPVStatus {
        input:
            coverage = SamtoolsCoverage.coverage
    }

    call GenotypeSNPsHuman {
        input:
            bam = SortAndIndexFinalBam.sorted_bam,
            bai = SortAndIndexFinalBam.sorted_bam_index,
            human_snp_targets_bed = human_snp_targets_bed,
            reference = reference,
            output_basename = output_basename
    }

    output {
        File final_bam = SortAndIndexFinalBam.sorted_bam
        File final_bam_index = SortAndIndexFinalBam.sorted_bam_index
        File umi_group_data = MergeBAMsAndGroupUMIs.umi_group_data
        File umi_duplication_metrics = CollectUMIDuplicationMetrics.umi_duplication_metrics
        File vcf = GenotypeSNPsHuman.vcf
        File coverage = SamtoolsCoverage.coverage
        String top_hpv_contig = DetermineHPVStatus.top_hpv_contig
        Int top_hpv_num_reads = DetermineHPVStatus.top_hpv_num_reads
        Float top_hpv_coverage = DetermineHPVStatus.top_hpv_coverage
        Boolean is_hpv_positive = DetermineHPVStatus.is_hpv_positive
        String secondary_hpv_types = DetermineHPVStatus.secondary_hpv_types
        File fastp_report_html = TrimAndFilter.fastp_report_html
        File fastp_report_json = TrimAndFilter.fastp_report_json
        File pre_trimmed_r1_fastqc_html = PreTrimmedFastQC.r1_fastqc_html
        File pre_trimmed_r2_fastqc_html = PreTrimmedFastQC.r2_fastqc_html
        File post_trimmed_r1_fastqc_html = PostTrimmedFastQC.r1_fastqc_html
        File post_trimmed_r2_fastqc_html = PostTrimmedFastQC.r2_fastqc_html
        File pre_consensus_alignment_summary_metrics = PreConsensusAlignmentSummaryMetrics.alignment_summary_metrics
        File pre_consensus_flagstat = PreConsensusFlagstat.flagstat
        File pre_consensus_insert_size_metrics = PreConsensusInsertSizeMetrics.insert_size_metrics
        File pre_consensus_insert_size_plot = PreConsensusInsertSizeMetrics.insert_size_plot
        File pre_consensus_ontarget_reads = PreConsensusCountOnTargetReads.ontarget_reads
        File pre_consensus_hs_metrics = PreConsensusHybridSelectionMetrics.hs_metrics
        File pre_consensus_per_base_coverage = PreConsensusHybridSelectionMetrics.per_base_coverage
        File post_consensus_alignment_summary_metrics = PostConsensusAlignmentSummaryMetrics.alignment_summary_metrics
        File post_consensus_flagstat = PostConsensusFlagstat.flagstat
        File post_consensus_insert_size_metrics = PostConsensusInsertSizeMetrics.insert_size_metrics
        File post_consensus_insert_size_plot = PostConsensusInsertSizeMetrics.insert_size_plot
        File post_consensus_ontarget_reads = PostConsensusCountOnTargetReads.ontarget_reads
        File post_consensus_hs_metrics = PostConsensusHybridSelectionMetrics.hs_metrics
        File post_consensus_per_base_coverage = PostConsensusHybridSelectionMetrics.per_base_coverage
    }
}