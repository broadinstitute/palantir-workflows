# HPV DeepSeek

This directory contains the WDL implementation of the HPV DeepSeek assay, a UMI/duplex-sequencing based
targeted panel for detecting and characterizing human papillomavirus (HPV) cell-free DNA in plasma. The
panel baits both HPV genomes and a set of human (hg38) targets, and is sequenced with duplex UMIs so that
very low allele-fraction HPV signal can be separated from sequencing error. The workflows here take raw
paired-end FASTQs through UMI-aware consensus alignment, HPV genotype calling, somatic variant calling on
the human targets, HPV-specific tertiary analyses (integration breakpoints, HPV16 sublineage assignment,
high-risk SNPs), and a normalization step that converts HPV depth into a per-mL-plasma quantity.

The reference used throughout is a combined reference containing both the hg38 human contigs and one contig
per HPV genotype (HPV contigs are identified by their names starting with `HPV`, and the HPV16 reference
contig is named `HPV16_Ref`). Building that reference, the bait/target interval lists, and the various
resource files is outside the scope of these WDLs — they are all supplied as inputs.

[HPVDeepSeek.wdl](HPVDeepSeek.wdl) is the entry point most users should run: it chains the four
subworkflows together and re-exports all of their outputs. The subworkflows can also be run standalone if
you already have the intermediate BAMs. All five WDLs are registered in Dockstore via
[/.dockstore.yml](../.dockstore.yml) under the names `HPVDeepSeek`, `HPVDeepSeekGenotyping`,
`HPVDeepSeekSomaticVariantCalling`, `HPVDeepSeekTertiaryAnalysis`, and `HPVDeepSeekNormalization`.

This directory contains the following WDLs:
- [HPVDeepSeek](#hpvdeepseek): Top-level workflow; runs genotyping, somatic variant calling, tertiary analysis, and normalization end-to-end from FASTQs
- [HPVDeepSeekGenotyping](#hpvdeepseekgenotyping): FASTQ QC/trimming, UMI extraction, alignment, simplex and duplex consensus calling, panel metrics, HPV genotype status, and human SNP genotyping
- [HPVDeepSeekSomaticVariantCalling](#hpvdeepseeksomaticvariantcalling): Mutect2 tumor-only somatic calling on the duplex consensus BAM, with contamination, orientation-bias, BLAST mapping and population-AF filters, plus Funcotator annotation
- [HPVDeepSeekTertiaryAnalysis](#hpvdeepseektertiaryanalysis): HPV integration breakpoint detection, HPV16 sublineage phylogenetic assignment, and high-risk HPV SNP reporting
- [HPVDeepSeekNormalization](#hpvdeepseeknormalization): Normalizes HPV depth against human depth and converts it to an HPV quantity per mL of plasma

### Docker images

These workflows use a mixture of public Broad images and private images in the Hydrogen Artifact Registry
(`us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/...`). All image references are
hardcoded in the task `runtime` blocks — none of them are exposed as workflow inputs, so changing a tool
version requires editing the WDL.

Public images:
- `us.gcr.io/broad-gatk/gatk:4.6.2.0` — the GATK tasks in somatic variant calling
- `us.gcr.io/broad-dsde-methods/bcftools:v1.4` — `bcftools` genotyping tasks
- `us.gcr.io/broad-dsde-methods/liquidbiopsy:0.0.3.5` — the BLAST-based mapping filter (also provides the
  default `filter_alt_ref_positions.py` script and `blastn` binary paths)
- `us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian` — UMI duplication metric computation
- `gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10` — FastQC

Hydrogen images (pinned by digest where applicable):
- `.../kockan/hds` — the general-purpose image carrying `gatk`, `fgbio`, `bwa`, `samtools`, and `fastp`
- `.../kockan/fgbio` — a separate fgbio image used for duplex consensus calling and duplex-seq metrics
- `.../kockan/gatk-custom:v4.6.2.0a` — custom GATK build used for Funcotator
- `.../kockan/simple_pysam` — python/pysam/pandas image for HPV status and normalization
- `.../kockan/breakpoint_detector` — carries `/breakpoint_detector-v3.7.py`
- `.../kockan/hpv_sublineage` — carries `bcftools`, `samtools`, `muscle`, `seqret` (EMBOSS), and `phyml`
- `.../kockan/hpv_toytree` — carries `toytree` and `toyplot`

### Testing

There are currently no automated test configurations for these workflows: `HPVDeepSeek` does not appear in
[/test/watt_config.yml](../test/watt_config.yml), and there are no example input JSONs checked into this
repository. Inputs must be assembled by hand (or from a Terra workspace configuration).

## HPVDeepSeek

### Summary

Top-level workflow and the recommended entry point. It runs, in order:

1. [HPVDeepSeekGenotyping](#hpvdeepseekgenotyping) on the input FASTQ pair, producing raw, simplex-consensus,
   and duplex-consensus BAMs plus panel/UMI metrics, `samtools coverage`, and the HPV status table.
2. [HPVDeepSeekSomaticVariantCalling](#hpvdeepseeksomaticvariantcalling) on the **duplex** consensus BAM,
   using `hg38_target_interval_list` as the Mutect2 intervals.
3. [HPVDeepSeekTertiaryAnalysis](#hpvdeepseektertiaryanalysis) on the **simplex** consensus BAM.
4. [HPVDeepSeekNormalization](#hpvdeepseeknormalization) on the simplex consensus BAM plus the HPV status
   table from step 1.

Steps 2–4 do not depend on each other, so Cromwell will run them concurrently once genotyping completes.
The `output_basename` input is used as the basename for essentially every output file, and is also passed
through as the `sample_id` for normalization.

### Inputs

Genotyping inputs:
- **String output_basename**: Basename used for all output files; also used as the sample id in the
  normalization output.
- **File r1_fastq**, **File r2_fastq**: Paired-end FASTQs (gzipped or not) for the sample.
- **File human_snp_targets_bed**: BED of human SNP positions to genotype with `bcftools mpileup`/`call`
  (used for fingerprinting/identity checks).
- **File reference**, **File reference_fai**, **File reference_dict**: Combined hg38 + HPV reference FASTA
  and its index and sequence dictionary.
- **File bwa_idx_amb**, **File bwa_idx_ann**, **File bwa_idx_bwt**, **File bwa_idx_pac**, **File bwa_idx_sa**:
  BWA index files for `reference`. These are localized next to the FASTA; they are not referenced directly
  in the command lines.
- **File hpv_bait_interval_list**, **File hpv_target_interval_list**: Picard-style interval lists over the HPV
  contigs, for `CollectHsMetrics`.
- **File hg38_bait_interval_list**, **File hg38_target_interval_list**: Picard-style interval lists over the
  human targets, for `CollectHsMetrics`. `hg38_target_interval_list` is also used as the Mutect2 interval
  list in somatic variant calling.
- **File low_risk_hpv_genotypes**: Plain text file with one HPV contig/genotype name per line. Detected
  genotypes in this list are flagged as not reportable in the HPV status table.
- **String bait_set_name**: Bait set name recorded in the `CollectHsMetrics` output.
- **String read_group_id**, **String read_group_sample_name**: Read group ID and sample name written into the
  unmapped BAM and the BWA `@RG` line.
- **String read_group_library_name = "LB_DEFAULT"**, **String read_group_platform = "ILLUMINA"**,
  **String read_group_platform_unit = "PU_DEFAULT"**, **String read_group_description = "DS_DEFAULT"**:
  Remaining read group fields. The placeholder defaults should be overridden for real production runs.
- **String read_structure = "3M2S+T"**: fgbio read structure applied to both reads when extracting UMIs —
  3 bases of UMI, 2 skipped bases, remainder is template.

Somatic variant calling inputs:
- **File gnomad**, **File gnomad_idx**: gnomAD germline resource VCF for Mutect2.
- **File pon**, **File pon_idx**: Panel of normals VCF for Mutect2.
- **File variants_for_contamination**, **File variants_for_contamination_idx**: Common-variant VCF used by
  `GetPileupSummaries`/`CalculateContamination`.
- **File realignment_index_bundle**: BWA-MEM index image used by `FilterAlignmentArtifacts`. Required even
  when `run_alignment_artifact_filter` is false.
- **String mapping_filter_python_script = "/usr/filter_alt_ref_positions.py"**: Path *inside the
  liquidbiopsy docker image* to the BLAST-based mapping filter script.
- **File blastdb_nhr**, **File blastdb_nin**, **File blastdb_nsq**: BLAST database files localized for the
  mapping filter task.
- **String blastn_path = "/usr/blastn_2.2.30+"**: Path inside the liquidbiopsy docker image to the `blastn`
  binary handed to the mapping filter script.
- **File funcotator_data_source**: Funcotator data sources tarball (`.tar.gz`).
- **Boolean run_alignment_artifact_filter = false**: If true, runs GATK `FilterAlignmentArtifacts` between
  `FilterMutectCalls` and the mapping filter.

Tertiary analysis inputs:
- **File high_risk_snps_hpv**: File whose first column is a list of HPV positions considered high risk; the
  workflow reports which of these are present in the sample's HPV16 variant calls.
- **File hpv16_sublineages**: Multi-FASTA of HPV16 sublineage reference sequences, concatenated with the
  sample consensus before multiple sequence alignment and tree building.

Normalization inputs:
- **File fp_intervals**: 4-column headerless BED (`chromosome`, `start`, `end`, `info`) of human intervals
  over which mean depth is computed to derive the human background depth.
- **Float ul_plasma**: Volume of plasma in microliters that the library was made from (converted internally
  to mL).
- **Float ng_cfdna**: Mass of cfDNA in nanograms that went into the library.

### Outputs

All outputs are pass-throughs of the corresponding subworkflow outputs; see the subworkflow sections for
descriptions.

- From [HPVDeepSeekGenotyping](#hpvdeepseekgenotyping): `raw_bam`, `raw_bam_index`, `simplex_bam`,
  `simplex_bam_index`, `duplex_bam`, `duplex_bam_index`, `simplex_umi_grouped_bam`, `simplex_umi_group_data`,
  `duplex_umi_grouped_bam`, `duplex_umi_group_data`, `simplex_umi_duplication_metrics`,
  `duplex_umi_duplication_metrics`, `vcf`, `coverage`, `hpv_status`, `fastp_report_html`,
  `fastp_report_json`, `pre_trimmed_r1_fastqc_html`, `pre_trimmed_r2_fastqc_html`,
  `post_trimmed_r1_fastqc_html`, `post_trimmed_r2_fastqc_html`,
  `pre_consensus_alignment_summary_metrics`, `pre_consensus_insert_size_metrics`,
  `pre_consensus_insert_size_plot`, `post_consensus_alignment_summary_metrics`,
  `post_consensus_insert_size_metrics`, `post_consensus_insert_size_plot`, `raw_hpv_hs_metrics`,
  `raw_hpv_per_target_coverage`, `raw_hg38_hs_metrics`, `raw_hg38_per_target_coverage`,
  `simplex_hpv_hs_metrics`, `simplex_hpv_per_target_coverage`, `simplex_hg38_hs_metrics`,
  `simplex_hg38_per_target_coverage`, `duplex_hpv_hs_metrics`, `duplex_hpv_per_target_coverage`,
  `duplex_hg38_hs_metrics`, `duplex_hg38_per_target_coverage`, `family_sizes`, `duplex_family_sizes`,
  `duplex_yield_metrics`, `umi_counts`, `duplex_qc`
- From [HPVDeepSeekSomaticVariantCalling](#hpvdeepseeksomaticvariantcalling): `contamination_table`,
  `unfiltered_vcf`, `unfiltered_vcf_idx`, `mutect2_stats`, `filter_mutect_calls_stats`, `filtered_vcf`,
  `filtered_vcf_idx`, `funcotated_maf`
- From [HPVDeepSeekTertiaryAnalysis](#hpvdeepseektertiaryanalysis): `analysis_log`, `breakpoints`,
  `detailed_integration_summary`, `integration_breakpoints`, `integration_summary`,
  `multiple_sequence_alignment`, `phylip_formatted_msa`, `phylogenetic_tree_stats`, `phylogenetic_tree`,
  `phylogenetic_tree_visualization`, `sublineage_call`, `high_risk_snps_found`
- From [HPVDeepSeekNormalization](#hpvdeepseeknormalization): `normalized_hpv`

## HPVDeepSeekGenotyping

### Summary

Takes raw paired-end FASTQs to duplex- and simplex-consensus BAMs, collects QC across the whole process, and
calls the sample's HPV genotype status.

Steps, in order:

1. **FastQC** on the input FASTQs (pre-trimming).
2. **`gatk FastqToSam`** to build an unmapped BAM with the supplied read group.
3. **`fgbio ExtractUmisFromBam`** using `read_structure` on both reads, writing the UMI to the `RX` tag and
   appending it to the read name.
4. **`gatk SamToFastq`** back to FASTQ, then **`fastp`** trimming/filtering (auto adapter detection, 5-base
   quality sliding window at Q20, drop reads with >40% low-quality bases, poly-G trimming, minimum length 75,
   adapter trimming when only one end matches). Produces the fastp HTML/JSON reports.
5. **FastQC** again on the trimmed FASTQs (post-trimming).
6. **`bwa mem`** (`-M`, `-K 100000000`) with the read group, piped through `samtools view`, then
   `samtools sort`/`index`. This sorted BAM is the `raw_bam` output.
7. Pre-consensus QC on `raw_bam`: `CollectAlignmentSummaryMetrics`, `CollectInsertSizeMetrics`, and
   `CollectHsMetrics` against both the HPV and hg38 bait/target interval lists.
8. **`gatk MergeBamAlignment`** of the aligned BAM with the UMI-extracted unmapped BAM (queryname sorted,
   aligned reads only, `MostDistant` primary alignment strategy).
9. **`samtools view -f 2 -q 1`** (properly paired, MAPQ ≥ 1) followed by **`fgbio GroupReadsByUmi`**, run
   twice: once with the `adjacency` strategy (simplex) and once with the `paired` strategy (duplex), both
   with `--edits 1` on the `RX` tag. Each produces a UMI-grouped BAM and a family size histogram.
10. UMI duplication metrics computed from each family size histogram (percent duplication and estimated
    library size, derived from total vs. unique fragments).
11. **`fgbio CallMolecularConsensusReads`** on the simplex grouping and **`fgbio CallDuplexConsensusReads`**
    on the duplex grouping (both `--min-reads 1`, max 50 reads (per strand for duplex),
    `--min-input-base-quality 20`, pre/post-UMI error rates 45/40).
12. **`fgbio CollectDuplexSeqMetrics`** on the duplex UMI-grouped BAM.
13. Both consensus BAMs are converted back to FASTQ, realigned with `bwa mem` (this time with `-Y` soft
    clipping of supplementary alignments), queryname sorted, merged back with their unmapped consensus BAMs
    (retaining `RX`, coordinate sorted), filtered again with `samtools view -f 2 -q 1`, and sorted/indexed.
    These are the `simplex_bam` and `duplex_bam` outputs.
14. Post-consensus QC on the simplex BAM (`CollectAlignmentSummaryMetrics`, `CollectInsertSizeMetrics`) and
    `CollectHsMetrics` on both consensus BAMs against both interval list pairs.
15. **`samtools coverage`** on the duplex BAM, then HPV status determination: contigs whose names start with
    `HPV` and that have at least 2 duplex reads are reported, with an `Is_Reportable` flag that is false for
    genotypes listed in `low_risk_hpv_genotypes`.
16. **`bcftools mpileup | bcftools call -mv`** on the simplex BAM restricted to `human_snp_targets_bed`,
    producing the human SNP genotype VCF.

### Inputs

Same as the genotyping inputs of [HPVDeepSeek](#hpvdeepseek), with two differences:
- **String read_structure** has no default here and must be supplied (the top-level workflow defaults it to
  `"3M2S+T"`).
- The read group defaults differ: **read_group_library_name = "LB_TEST"**,
  **read_group_platform_unit = "PU_TEST"**, **read_group_description = "KAPA_TE"**.

### Outputs

BAMs:
- **File raw_bam** / **File raw_bam_index**: Pre-consensus coordinate-sorted alignment of the trimmed reads.
- **File simplex_bam** / **File simplex_bam_index**: Simplex (single-strand) consensus BAM, realigned,
  filtered, coordinate sorted and indexed.
- **File duplex_bam** / **File duplex_bam_index**: Duplex consensus BAM, same treatment.
- **File simplex_umi_grouped_bam**, **File duplex_umi_grouped_bam**: `GroupReadsByUmi` output BAMs.

Metrics and QC:
- **File simplex_umi_group_data**, **File duplex_umi_group_data**: fgbio family size histograms.
- **File simplex_umi_duplication_metrics**, **File duplex_umi_duplication_metrics**: Two-row TSVs with
  `PERCENT_DUPLICATION` and `ESTIMATED_LIBRARY_SIZE`.
- **File fastp_report_html**, **File fastp_report_json**: fastp trimming reports.
- **File pre_trimmed_r1_fastqc_html**, **File pre_trimmed_r2_fastqc_html**,
  **File post_trimmed_r1_fastqc_html**, **File post_trimmed_r2_fastqc_html**: FastQC reports.
- **File pre_consensus_alignment_summary_metrics**, **File pre_consensus_insert_size_metrics**,
  **File pre_consensus_insert_size_plot**: Picard metrics on `raw_bam`.
- **File post_consensus_alignment_summary_metrics**, **File post_consensus_insert_size_metrics**,
  **File post_consensus_insert_size_plot**: Picard metrics on `simplex_bam`.
- **File raw_hpv_hs_metrics** / **raw_hpv_per_target_coverage**, **raw_hg38_hs_metrics** /
  **raw_hg38_per_target_coverage**, **simplex_hpv_hs_metrics** / **simplex_hpv_per_target_coverage**,
  **simplex_hg38_hs_metrics** / **simplex_hg38_per_target_coverage**, **duplex_hpv_hs_metrics** /
  **duplex_hpv_per_target_coverage**, **duplex_hg38_hs_metrics** / **duplex_hg38_per_target_coverage**:
  `CollectHsMetrics` output for each combination of BAM (raw/simplex/duplex) and bait set (HPV/hg38).
  Coverage cap is 100000.
- **File family_sizes**, **File duplex_family_sizes**, **File duplex_yield_metrics**, **File umi_counts**,
  **File duplex_qc**: `fgbio CollectDuplexSeqMetrics` output (the last is a PDF).

Calls:
- **File coverage**: `samtools coverage` output on the duplex BAM, all contigs.
- **File hpv_status**: TSV with columns `HPV_Genotype`, `Num_Duplex_Reads`, `%_Genomic_Coverage`,
  `Is_Reportable`, one row per detected HPV contig.
- **File vcf**: Human SNP genotype VCF from `bcftools call` over `human_snp_targets_bed`.

### Notes

`CollectHsMetrics` is run with a fixed 512 GB SSD and most upstream tasks default to SSD disks of 512 GB or
1 TB regardless of input size, so this workflow is disk-expensive relative to the amount of data a targeted
panel produces. The `use_ssd`/`min_ssd_size_gb`/`disk_size_gb` task-level inputs are not exposed at the
workflow level, so tuning them requires editing the WDL or using task-level input overrides.

## HPVDeepSeekSomaticVariantCalling

### Summary

Tumor-only somatic variant calling with Mutect2 on the duplex consensus BAM, restricted to the human target
intervals. Steps, in order:

1. **`gatk CollectSequencingArtifactMetrics`** on the input BAM. Note this task's output is not consumed by
   any downstream task and is not a workflow output — it is marked `# !UnusedCall` in the source.
2. **`gatk Mutect2`** in tumor-only mode with the germline resource and panel of normals, over
   `mutect_target_intervals`. Notable settings tuned for a deep, duplex-consensus panel:
   `--read-filter NotSupplementaryAlignmentReadFilter`, `--dont-use-soft-clipped-bases true`,
   `--af-of-alleles-not-in-resource 0.001`, `--tumor-lod-to-emit 0`, `--initial-tumor-lod 0`,
   `--max-reads-per-alignment-start 0` (no downsampling), and `--pcr-snv-qual 70 --pcr-indel-qual 70`.
   The same task then runs **`gatk GetPileupSummaries`** against `variants_for_contamination`; the task
   preserves and exits with the Mutect2 exit code so a pileup failure does not mask a Mutect2 failure.
3. **`gatk LearnReadOrientationModel`** on the f1r2 counts and **`gatk CalculateContamination`** on the
   pileup summaries (also producing tumor segmentation).
4. **`gatk FilterMutectCalls`** using the contamination table, tumor segmentation, orientation-bias priors,
   and Mutect2 stats.
5. Optionally (`run_alignment_artifact_filter`), **`gatk FilterAlignmentArtifacts`** using the BWA-MEM index
   image. When enabled, its output becomes the input to the following steps.
6. **Mapping filter**: `filter_alt_ref_positions.py` (python2.7, from the liquidbiopsy image) BLASTs the
   variant-supporting positions against the provided BLAST database and emits a VCF of positions that pass.
   The result is bgzipped and tabix indexed.
7. **`gatk VariantFiltration`**, run twice: first flagging `POPAF < 3.0` (phred-scaled, i.e. population
   allele frequency > 0.001) as `germline`, then applying the mapping filter VCF as a mask with
   `--filter-not-in-mask true` and mask name `mapping_filter`, so variants *not* found in the mapping filter
   output get flagged.
8. **`gatk Funcotator`** on the filtered VCF, `--ref-version hg38`, MAF output, `BEST_EFFECT` transcript
   selection, `--prefer-mane-transcripts`, `--remove-filtered-variants` (filtered variants dropped), gnomAD
   data sources disabled.

### Inputs

- **String output_basename**: Basename for all output files.
- **File tumor_bam**, **File tumor_bai**: Input BAM and index. The top-level workflow passes the **duplex**
  consensus BAM here.
- **File mutect_target_intervals**: Intervals for Mutect2 and `GetPileupSummaries`. The top-level workflow
  passes `hg38_target_interval_list`.
- **File reference**, **File reference_fai**, **File reference_dict**: Combined hg38 + HPV reference.
- **File gnomad**, **File gnomad_idx**: Germline resource for Mutect2.
- **File pon**, **File pon_idx**: Panel of normals for Mutect2.
- **File variants_for_contamination**, **File variants_for_contamination_idx**: Sites used by
  `GetPileupSummaries`.
- **File realignment_index_bundle**: BWA-MEM index image for `FilterAlignmentArtifacts`. Required as an
  input even when the filter is disabled.
- **String mapping_filter_python_script = "/usr/filter_alt_ref_positions.py"**: In-container path to the
  mapping filter script.
- **File blastdb_nhr**, **File blastdb_nin**, **File blastdb_nsq**: BLAST database files. They are localized
  so `blastn` can find them, but are not named on the command line.
- **String blastn_path = "/usr/blastn_2.2.30+"**: In-container path to the `blastn` binary.
- **File funcotator_data_source**: Funcotator data sources `.tar.gz`; extracted in the task.
- **Boolean run_alignment_artifact_filter = false**: Toggle for `FilterAlignmentArtifacts`.

### Outputs

- **File contamination_table**: `CalculateContamination` output.
- **File unfiltered_vcf** / **File unfiltered_vcf_idx**: Raw Mutect2 calls.
- **File mutect2_stats**: Mutect2 `.stats` file.
- **File filter_mutect_calls_stats**: `FilterMutectCalls` filtering statistics.
- **File filtered_vcf** / **File filtered_vcf_idx**: Final VCF after `FilterMutectCalls` (and optionally
  `FilterAlignmentArtifacts`), plus the `germline` POPAF filter and the `mapping_filter` mask filter.
  Filters are annotated, not removed.
- **File funcotated_maf**: Funcotator MAF (`<output_basename>.maf`) with filtered variants removed.

### Notes

`CalculateContamination` and the germline POPAF filter both assume the panel includes enough human common
sites to be informative; with an HPV-heavy library the contamination estimate may be based on very few
sites. The `Funcotate` task exposes many optional inputs (transcript lists, annotation defaults/overrides,
excluded fields, interval list, extra args) at the task level, but the workflow hardcodes the values listed
above and does not surface them as workflow inputs.

## HPVDeepSeekTertiaryAnalysis

### Summary

HPV-specific downstream analyses. All four analyses run on the same input BAM (the top-level workflow passes
the **simplex** consensus BAM) and are independent of each other except for the tree drawing step.

1. **`DetectHPVIntegrationBreakpoints`** runs `/breakpoint_detector-v3.7.py` from the `breakpoint_detector`
   image, which scans the BAM for HPV–human junctions and writes an analysis log, a breakpoints table, an
   integration breakpoints table, and summary/detailed-summary tables.
2. **`Sublineages`** counts reads on `HPV16_Ref` with `samtools view -c`. If there are fewer than
   `read_threshold` (default 10) reads, the sample is called `HPV16_negative` with a `NA` distance and the
   alignment/tree outputs are created empty. Otherwise it builds a consensus HPV16 sequence
   (`bcftools mpileup --max-depth 8000 | bcftools call -mv | bcftools norm`, then
   `samtools faidx ... | bcftools consensus`), concatenates it with the `hpv16_sublineages` reference
   FASTA, aligns with `muscle`, converts to PHYLIP with EMBOSS `seqret`, and builds a tree with `phyml`.
3. **`SublineagesDrawTree`** loads the phyml tree with `toytree`, computes the tip-to-tip patristic distance
   matrix, and reports the sublineage closest to `HPV16_Ref` along with the distance, plus a PDF rendering
   of the tree. If the tree file is empty (the HPV16-negative case), it passes the upstream sublineage call
   through unchanged and emits an empty PDF.
4. **`HPVHighRiskSNPs`** runs `bcftools mpileup --max-depth 800000 | bcftools call -mv` restricted to
   `HPV16_Ref`, builds a variant table with an allele fraction computed from the `AD` field, and reports
   which of the positions listed in `high_risk_snps_hpv` were found.

### Inputs

- **String output_basename**: Basename for all output files; also used as the `run_id`/`library_id` in the
  sublineage call CSV.
- **File tumor_bam**, **File tumor_bai**: Input BAM and index. The top-level workflow passes the **simplex**
  consensus BAM.
- **File high_risk_snps_hpv**: File whose first whitespace/tab-delimited column lists high-risk HPV
  positions to look for in the variant table.
- **File reference**: Combined hg38 + HPV reference FASTA, used for the HPV16 pileups and consensus.
- **File hpv16_sublineages**: Multi-FASTA of HPV16 sublineage reference sequences.

The `Sublineages` task additionally takes **Int read_threshold = 10** (minimum `HPV16_Ref` read count before
attempting sublineage assignment) and `HPVHighRiskSNPs` takes **String genotype = "HPV16_Ref"** (the contig
to pile up), but neither is exposed at the workflow level.

### Outputs

Integration:
- **File analysis_log**: Breakpoint detector log.
- **File breakpoints**: All detected breakpoints.
- **File integration_breakpoints**: Breakpoints classified as HPV integration events.
- **File integration_summary**, **File detailed_integration_summary**: Summary tables of the integration
  calls.

Sublineage:
- **File multiple_sequence_alignment**: MUSCLE alignment of the sample HPV16 consensus with the sublineage
  references (`.afa`). Empty if HPV16-negative.
- **File phylip_formatted_msa**: The same alignment in PHYLIP format. Empty if HPV16-negative.
- **File phylogenetic_tree_stats**, **File phylogenetic_tree**: phyml stats and Newick tree. Empty if
  HPV16-negative.
- **File phylogenetic_tree_visualization**: PDF rendering of the tree. Empty if HPV16-negative.
- **File sublineage_call**: CSV with header `run_id,library_id,closest_sublineage,patristic_distance`. For
  HPV16-negative samples `closest_sublineage` is `HPV16_negative` and the distance is `NA`.

High-risk SNPs:
- **File high_risk_snps_found**: Tab-delimited list of high-risk positions found in the sample.

### Notes

Sublineage assignment is HPV16-only — the `HPV16_Ref` contig name is hardcoded in the task command lines.
The empty-file pattern in the HPV16-negative branch means downstream consumers must handle zero-byte
alignment/tree outputs rather than assuming they are always populated.

## HPVDeepSeekNormalization

### Summary

Converts raw HPV depth into a normalized quantity per mL of plasma, so that HPV load is comparable across
samples with different input amounts and different sequencing depths. Single task (`NormalizeHPV`), all in
python with `pysam` and `pandas`:

1. Read `fp_intervals` (a headerless 4-column BED) and append one interval per HPV contig found in the
   simplex BAM header, spanning the full contig.
2. For every interval, pile up the simplex BAM and count read bases whose `cD` tag (consensus depth) is at
   least 5, then divide by the number of pileup columns to get a mean depth. `max_depth` is set to 1,000,000
   and overlapping mate bases are ignored.
3. Compute the human background depth as the median of the interval mean depths over non-HPV, non-`chrX`,
   non-`chrY` contigs.
4. Keep only the HPV genotypes listed in the `hpv_status` table (i.e. the ones the genotyping workflow
   actually detected) and compute
   `HPV_Mean_Depth_Over_hg38_Median_Depth = mean_depth / hg38_median_depth` and
   `HPV_Quantity = HPV_Mean_Depth_Over_hg38_Median_Depth * ((ng_cfDNA / 0.0033) / mL_Plasma)`, where
   `0.0033` converts nanograms of cfDNA to genome equivalents (ng per haploid genome) and `mL_Plasma` is
   `ul_plasma / 1000`.

### Inputs

- **String sample_id**: Used as the output file basename. The top-level workflow passes `output_basename`.
- **File simplex_bam**, **File simplex_bam_index**: Simplex consensus BAM and index. The `cD` tag written by
  fgbio consensus calling must be present, so this must be a consensus BAM.
- **File hpv_status**: The `hpv_status` TSV from [HPVDeepSeekGenotyping](#hpvdeepseekgenotyping); only its
  `HPV_Genotype` column is used, to select which HPV contigs to report.
- **File fp_intervals**: Headerless 4-column BED (`chromosome`, `start`, `end`, `info`) of human intervals
  used to estimate the human background depth.
- **Float ul_plasma**: Plasma volume in microliters used for the library.
- **Float ng_cfdna**: cfDNA mass in nanograms used for the library.

### Outputs

- **File normalized_hpv**: `<sample_id>.normalized_hpv.tsv`, one row per detected HPV genotype with columns
  `HPV_Genotype`, `HPV_Mean_Depth_Over_hg38_Median_Depth`, `ng_cfDNA`, `mL_Plasma`, `HPV_Quantity`.

### Notes

The pileup loop is a pure-python per-base scan over every interval in `fp_intervals` plus every HPV contig,
so runtime scales with the total interval size; the task defaults to 2 CPUs and 16 GB. If no positions are
covered for an interval the mean depth is recorded as 0, which will make `HPV_Quantity` zero (or, if the
human background median is 0, undefined).
