version 1.0

# TODO: Add reference and bam indices + auxiliary files wherever required
# TODO: For GATK and fgbio, stick to consistent command-line argument convention after picking one
# TODO: Decide on input/output/intermediary file naming convensions based on the proposed structure of Sample Names

task VerifyPipelineInputs {
    input {
        File? bam
        File? r1_fastq
        File? r2_fastq

        Int? cpu = 1
        Int? memory_gb = 8
        Int? disk_size_gb = ceil(size(bam, "GiB") + size(r1_fastq,"GiB") + size(r2_fastq, "GiB")) + 32
    }

    command <<<
        set -e
        python3 <<CODE

        fastq_flag = 0
        bam = "~{bam}"
        r1_fastq = "~{r1_fastq}"
        r2_fastq = "~{r2_fastq}"

        if bam and not r1_fastq and not r2_fastq:
            pass
        elif r1_fastq and r2_fastq and not bam:
            fastq_flag += 1
        else:
            raise ValueError("Invalid Input. Input must be either ubam or a pair of fastqs")

        with open("output.txt", "w") as f:
            if fastq_flag == 1:
                f.write("true")
            # Remaining case is that only bam is defined
            else:
                f.write("false")
        CODE
    >>>

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    }

    output {
        Boolean fastq_run = read_boolean("output.txt")
    }
}

task FastQC {
    input {
        File r1_fastq
        File r2_fastq
        Int? cpu = 2
        Int? num_threads = 4
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2 * (size(r1_fastq, "GiB") + size(r2_fastq, "GiB"))) + 50)
    }

    String r1_fastq_name = sub(sub(basename(r1_fastq), "\\.fastq.gz$", ""), "\\.fq.gz$", "")
    String r2_fastq_name = sub(sub(basename(r2_fastq), "\\.fastq.gz$", ""), "\\.fq.gz$", "")

    command <<<
        fastqc ~{r1_fastq} ~{r2_fastq} \
        --threads ~{num_threads} \
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
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
    }
}

task FastqToUbam {
    input {
        File r1_fastq
        File r2_fastq
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * (size(r1_fastq, "GiB") + size(r2_fastq, "GiB"))) + 50)
    }

    String prefix = basename(r1_fastq, ".fastq.gz")

    command <<<
        gatk FastqToSam \
        -F1 ~{r1_fastq} \
        -F2 ~{r2_fastq} \
        -O ~{prefix}.unmapped.bam \
        -SM SAMPLE
    >>>

    output {
        File ubam = "~{prefix}.unmapped.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# Extract the UMI sequence from the first 3 bases of each read, skip the next 3 bases, append a 'T',
# and add the resulting UMI to the RX tag and read name in the output BAM.
#
# --read-group-tag RX: Tells fgbio to place the extracted UMI into the RX tag of each read.
# -r 3M3S+T 3M3S+T: Specifies the regular expressions used to extract the UMI from each read sequence:
# 3M3S+T means:
# 3M: Match the first 3 bases (UMI).
# 3S: Skip (soft clip) the next 3 bases.
# +T: Append this to the UMI (the sequence of "T").
# This regex is applied to both R1 and R2 (paired-end reads).
# -a true: Indicates that the UMI should be appended to the read name (QNAME) in addition to the RX tag.
task ExtractUMIs {
    input {
        File input_ubam
        String read_group_tag = "RX"
        String read_structure = "3M3S+T"
        String append_umi_to_qname = "true"
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(input_ubam, "GiB")) + 50)
    }

    String prefix = basename(input_ubam, ".bam")

    command <<<
        fgbio ExtractUmisFromBam \
        -i ~{input_ubam} \
        -o ~{prefix}.umi_extracted.bam \
        -r ~{read_structure} ~{read_structure} \
        -t ~{read_group_tag} \
        -a ~{append_umi_to_qname}
    >>>

    output {
        File umi_extracted_bam = "~{prefix}.umi_extracted.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# Convert UMI extracted uBAM to FASTQ before adapter trimming and filtering
task UmiExtractedBamToFastq {
    input {
        File umi_extracted_bam
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(umi_extracted_bam, "GiB")) + 50)
    }

    String prefix = basename(umi_extracted_bam, ".bam")

    command <<<
        gatk SamToFastq \
        --INPUT ~{umi_extracted_bam} \
        --FASTQ ~{prefix}_R1.fastq \
        --SECOND_END_FASTQ ~{prefix}_R2.fastq
    >>>

    output {
        File umi_extracted_fastq1 = "~{prefix}_R1.fastq"
        File umi_extracted_fastq2 = "~{prefix}_R2.fastq"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
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
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * (size(fastq1, "GiB") + size(fastq2, "GiB"))) + 50)
    }

    String prefix = basename(fastq1, ".umi_extracted_R1.fastq")

    command <<<
        fastp \
        -i ~{fastq1} \
        -I ~{fastq2} \
        -o ~{prefix}.trimmed_R1.fastq \
        -O ~{prefix}.trimmed_R2.fastq \
        -g \
        -W 5 \
        -q 20 \
        -u 40 \
        -3 \
        -l 75 \
        -c \
        -h ~{prefix}.fastp_report.html \
        -j ~{prefix}.fastp_report.json
    >>>

    output {
        File fastq1_trimmed = "~{prefix}.trimmed_R1.fastq"
        File fastq2_trimmed = "~{prefix}.trimmed_R2.fastq"
        File fastp_report_html = "~{prefix}.fastp_report.html"
        File fastp_report_json = "~{prefix}.fastp_report.json"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
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
        File fastq1
        File fastq2
        File reference
        File bwa_idx_amb
        File bwa_idx_ann
        File bwa_idx_bwt
        File bwa_idx_pac
        File bwa_idx_sa
        String read_group_id
        String read_group_sample
        String platform
        Int? cpu = 16
        Int? num_threads = 16
        Int? memory_gb = 64
        Int? disk_size_gb = ceil((3 * (size(fastq1, "GiB") + size(fastq2, "GiB"))) + size(reference, "GiB") + 100)
    }

    String prefix = basename(fastq1, ".unmapped.trimmed_R1.fastq")

    command <<<
        bwa mem \
        -t ~{num_threads} \
        -K 100000000 \
        -R '@RG\tID:~{read_group_id}\tDS:KAPA_TE\tPL:~{platform}\tLB:lib1\tSM:~{read_group_sample}\tPU:unit1' \
        -M \
        ~{reference} \
        ~{fastq1} \
        ~{fastq2} \
        | samtools view --threads ~{num_threads} -o ~{prefix}.bam -
    >>>

    output {
        File bam = "~{prefix}.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# Sort and index aligned BAM
task SortAndIndexBam {
    input {
        File bam
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
    }

    String prefix = basename(bam, ".bam")

    command <<<
        samtools sort \
        -o ~{prefix}.sorted.bam \
        -O bam \
        -T $~{prefix}.bam.temp \
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
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# GATK Sort BAM by queryname
task GATKSortBam {
    input {
        File bam
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
    }

    String prefix = basename(bam, ".bam")

    command <<<
        gatk SortSam \
        -I ~{bam} \
        -O ~{prefix}.sorted.bam \
        -SO queryname
    >>>

    output {
        File sorted_bam = "~{prefix}.sorted.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# 1. Merge BWA-aligned reads with original UMI-tagged reads to produce a BAM file that has alignment info + UMI tags, necessary for UMI-aware deduplication.
# 2. Group together reads that share the same UMI sequence (or a similar one within a defined edit distance) and originate from the same genomic location—used for UMI-aware deduplication.
# fgbio:
#--input: Input BAM file containing reads aligned to the reference with UMI tags (RX) and filtered.
#--output: Output BAM where each read is tagged with a group ID (MI tag) indicating its molecular family.
#--strategy=adjacency: Uses an adjacency graph to group UMIs—UMIs within one base mismatch (--edits=1) are grouped if one is more abundant.
#--edits=1: Allows grouping of UMIs within 1 base mismatch (to account for sequencing errors).
#-t RX: Specifies the tag where UMIs are stored (here, the standard RX tag).
#-f ...umi_group_data.txt: Outputs grouping statistics and metadata for each molecular family.
task MergeBAMsAndGroupUMIs {
    input {
        File aligned_bam
        File unmapped_umi_extracted_bam
        File reference
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(aligned_bam, "GiB") + size(unmapped_umi_extracted_bam, "GiB")) + 100)
    }

    String prefix = basename(aligned_bam, ".bam")

    command <<<
        gatk MergeBamAlignment \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ATTRIBUTES_TO_REMOVE NM \
        --ATTRIBUTES_TO_REMOVE MD \
        --ALIGNED_BAM ~{aligned_bam} \
        --UNMAPPED_BAM ~{unmapped_umi_extracted_bam} \
        --OUTPUT ~{prefix}.umi_extracted.aligned.merged.bam \
        --REFERENCE_SEQUENCE ~{reference} \
        --SORT_ORDER queryname \
        --ALIGNED_READS_ONLY true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --CLIP_OVERLAPPING_READS false

        fgbio GroupReadsByUmi \
        --input=~{prefix}.umi_extracted.aligned.merged.bam \
        --output=~{prefix}.umi_grouped.bam \
        --strategy=adjacency \
        --edits=1 \
        -t RX \
        -f ~{prefix}.umi_group_data.txt
    >>>

    output {
        File umi_grouped_bam = "~{prefix}.umi_grouped.bam"
        File umi_group_data = "~{prefix}.umi_group_data.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
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
        File consensus_aligned_bam
        File consensus_unmapped_bam
        File reference
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(consensus_aligned_bam, "GiB") + size(consensus_unmapped_bam, "GiB")) + 100)
    }

    String prefix = basename(consensus_aligned_bam, ".bam")

    command <<<
        gatk MergeBamAlignment \
        --ALIGNED_BAM ~{consensus_aligned_bam} \
        --UNMAPPED_BAM ~{consensus_unmapped_bam} \
        --OUTPUT ~{prefix}.deduped.bam \
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
        File deduped_bam = "~{prefix}.deduped.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
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
        File umi_grouped_bam
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(umi_grouped_bam, "GiB")) + 100)
    }

    String prefix = basename(umi_grouped_bam, ".bam")

    command <<<
        fgbio CallMolecularConsensusReads \
        --input ~{umi_grouped_bam} \
        --output ~{prefix}.umi_consensus.unmapped.bam \
        --error-rate-post-umi 40 \
        --error-rate-pre-umi 45 \
        --output-per-base-tags false \
        --min-reads 1 \
        --max-reads 50 \
        --min-input-base-quality 20 \
        --read-name-prefix='consensus'
    >>>

    output {
        File umi_consensus_unmapped_bam = "~{prefix}.umi_consensus.unmapped.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# Convert consensus BAM to FASTQ:
task ConsensusBamToFastq {
    input {
        File umi_consensus_unmapped_bam
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(umi_consensus_unmapped_bam, "GiB")) + 50)
    }

    String prefix = basename(umi_consensus_unmapped_bam, ".bam")

    command <<<
        gatk SamToFastq \
        INPUT=~{umi_consensus_unmapped_bam} \
        FASTQ=~{prefix}.consensus_unmapped_R1.fastq \
        SECOND_END_FASTQ=~{prefix}.consensus_unmapped_R2.fastq
    >>>

    output {
        File consensus_unmapped_fastq1 = "~{prefix}.consensus_unmapped_R1.fastq"
        File consensus_unmapped_fastq2 = "~{prefix}.consensus_unmapped_R2.fastq"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# Reads and Coverage Calculation and results generation
# samtools coverage deduped_sorted.bam | awk '{if(NR==1){printf "%s\t%s\t%s\n",$1,$4,$6} else {printf "%s\t%s\t%s\n",$1,$4,$6}}' > sample.coverage.txt
task SamtoolsCoverage {
    input {
        File bam
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2 * size(bam, "GiB")) + 50)
    }

    String prefix = basename(bam, ".bam")

    command <<<
        samtools coverage ~{bam} | awk '{if(NR==1){printf "%s\t%s\t%s\n",$1,$4,$6} else {printf "%s\t%s\t%s\n",$1,$4,$6}}' > ~{prefix}.coverage.txt
    >>>

    output {
        File coverage = "~{prefix}.coverage.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

# TODO: The document doesn't have an implementation for the classification flowchart below. Add it ourselves
# HPV+ Classification
# Sample is considered HPV+ if both thresholds are met:
# Read count ≥ 10
# Percentage of HPV genome covered by read alignment ≥ 10%

# Human SNP Genotyping
task GenotypeSNPsHuman {
    input {
        File bam
        File human_snp_targets_bed
        File reference
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((2.5 * size(bam, "GiB")) + 50)
    }

    String prefix = basename(bam, ".bam")

    command <<<
        bcftools mpileup -f ~{reference} -R ~{human_snp_targets_bed} ~{bam} 2> ~{prefix}.mpileup.log | \
        bcftools call -mv -Ov -o ~{prefix}.vcf 2> ~{prefix}.call.log
    >>>

    output {
        File vcf = "~{prefix}.vcf"
        File mpileup_log = "~{prefix}.mpileup.log"
        File call_log = "~{prefix}.call.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

workflow HPVDeepSeek {
    input {
        File? ubam
        File? r1_fastq
        File? r2_fastq
        File human_snp_targets_bed
        File reference
        File reference_index
        File bwa_idx_amb
        File bwa_idx_ann
        File bwa_idx_bwt
        File bwa_idx_pac
        File bwa_idx_sa
        String read_group_id
        String read_group_sample
        String platform = "ILLUMINA"
    }

    call VerifyPipelineInputs {
        input:
            bam = ubam,
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq
    }

    if(VerifyPipelineInputs.fastq_run) {
        call FastQC {
            input:
                r1_fastq = select_first([r1_fastq]),
                r2_fastq = select_first([r2_fastq])
        }

        call FastqToUbam {
            input:
                r1_fastq = select_first([r1_fastq]),
                r2_fastq = select_first([r2_fastq])
        }
    }

    File ubam_to_use = select_first([ubam, FastqToUbam.ubam])

    call ExtractUMIs {
         input:
            input_ubam = ubam_to_use
    }

    call UmiExtractedBamToFastq {
        input:
            umi_extracted_bam = ExtractUMIs.umi_extracted_bam
    }

    call TrimAndFilter {
        input:
            fastq1 = UmiExtractedBamToFastq.umi_extracted_fastq1,
            fastq2 = UmiExtractedBamToFastq.umi_extracted_fastq2
    }

    call BwaMem {
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
            read_group_sample = read_group_sample,
            platform = platform
    }

    call SortAndIndexBam {
         input:
            bam = BwaMem.bam
    }

    call MergeBAMsAndGroupUMIs {
         input:
            aligned_bam = SortAndIndexBam.sorted_bam,
            unmapped_umi_extracted_bam = ExtractUMIs.umi_extracted_bam,
            reference = reference
    }

    call CallMolecularConsensusReads {
        input:
            umi_grouped_bam = MergeBAMsAndGroupUMIs.umi_grouped_bam
    }

    call ConsensusBamToFastq {
        input:
            umi_consensus_unmapped_bam = CallMolecularConsensusReads.umi_consensus_unmapped_bam
    }

    # Align consensus reads to the reference genome:
    # bwa mem -t THREADS \
    # -R '@RG\tID:A\tDS:KAPA_TE\tPL:ILLUMINA\tLB:lib1\tSM:sample1\tPU:unit1' \
    # -v 3 -Y -M -K 100000000 \
    # reference.fasta \
    # consensus_unmapped_R1.fastq consensus_unmapped_R2.fastq | \
    # samtools view -bh - > consensus_mapped_unsorted.bam
    call BwaMem as BwaMemRound2{
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
            read_group_sample = read_group_sample,
            platform = platform
    }

    call GATKSortBam as GATKSortBamConsensusAligned {
        input:
            bam = BwaMemRound2.bam
    }

    call GATKSortBam as GATKSortBamConsensusUnmapped{
        input:
            bam = CallMolecularConsensusReads.umi_consensus_unmapped_bam
    }

    call MergeConsensus {
         input:
            consensus_aligned_bam = GATKSortBamConsensusAligned.sorted_bam,
            consensus_unmapped_bam = GATKSortBamConsensusUnmapped.sorted_bam,
            reference = reference
    }

    call SortAndIndexBam as SortAndIndexFinalBam {
        input:
            bam = MergeConsensus.deduped_bam
    }

    call SamtoolsCoverage {
         input:
            bam = SortAndIndexFinalBam.sorted_bam
    }

    call GenotypeSNPsHuman {
        input:
            bam = SortAndIndexFinalBam.sorted_bam,
            human_snp_targets_bed = human_snp_targets_bed,
            reference = reference
    }

    output {
        File final_bam = SortAndIndexFinalBam.sorted_bam
        File final_bam_index = SortAndIndexFinalBam.sorted_bam_index
        File umi_group_data = MergeBAMsAndGroupUMIs.umi_group_data
        File vcf = GenotypeSNPsHuman.vcf
        File mpileup_log = GenotypeSNPsHuman.mpileup_log
        File call_log = GenotypeSNPsHuman.call_log
        File coverage = SamtoolsCoverage.coverage
        File fastp_report_html = TrimAndFilter.fastp_report_html
        File fastp_report_json = TrimAndFilter.fastp_report_json
        File? r1_fastqc_html = FastQC.r1_fastqc_html
        File? r2_fastqc_html = FastQC.r2_fastqc_html
    }
}