# Benchmark Phasing

This directory contains some workflows related to phasing. The concept of phasing includes understanding how nearby 
heterozygous variants are oriented relative to each other, i.e. on the same contig or opposite in diploid samples.
Here we have a workflow for computing this information, and for benchmarking samples that have already been phased.
When used together, this creates a pipeline for testing how phasing information has changed from one set of experimental
conditions to another.

1. [The Phasing WDL](#the-phasing-wdl)
2. [The Benchmarking WDL](#the-benchmarking-wdl)


## The Phasing WDL

This WDL adds phasing information to a VCF using read information.

The WDL [PhaseVCF](PhaseVCF.wdl) has the following input/output schema.

### Inputs

- `input_vcf`: the VCF to add phasing information to
- `input_vcf_index`: the index for the input VCF
- `input_sample_name`: optional but recommended input 
- `input_bam`: the read alignments (BAM) to use for linking nearby variants
- `input_bam_index`: the index for the input BAM
- `reference_fasta`: reference file 
- `reference_fasta_index`: index for the reference file

### Outputs

- `phased_vcf`: the input VCF with phasing information added, when possible
- `phased_vcf_index`: the index for the phased VCF
- `phase_stats`: a tsv of stats produced via `whatshap stats`

### Tool Options

Some options can be toggled in the `WhatsHapPhase` task for different behaviors. See the 
[tool documentation](https://whatshap.readthedocs.io/en/latest/guide.html) for a full description. Included options are:

- `ignore_read_groups`: (default = `true`) If false, phase samples using RGs
- `distrust_genotypes`: (default = `false`) If true, allow hets to be treated as homs if it allows for an optimal phasing solution
- `pedigree`: (optional) When included, use `.ped` file to strengthen phasing capabilities by leveraging trios consensus
- `tag`: (default = `PS`) tag for labeling phasing information; other option is `HP`
- `recombrate`: (optional) Toggle for non-human samples to better reflect species recombination rate


## The Benchmarking WDL

This WDL compares two phased VCFs and collects statistics about errors in the "callset" VCF against a "baseline" VCF. 
In contrast with benchmarking short variants for their genotypes, different error rate statistics are collected based
on inconsistencies in the call VCF's phasing against the baseline. More details can be found in the 
[WhatsHap documentation](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-compare).

The WDL [BenchmarkPhasing](BenchmarkPhasing.wdl) has the following input/output schema.

### Inputs

- `base_vcf`: the VCF to use as a baseline or "truth" data for phasing
- `base_vcf_index`: the index for the baseline VCF
- `base_sample_name`: (optional) the nickname to use for the sample in the `base_vcf`, as referenced in output tables
- `call_vcf`: the VCF "callset" to compare phasing against the baseline
- `call_vcf_index`: the index for the call VCF
- `call_sample_name`: (optional) the nickname to use for the sample in the `call_vcf`, as referenced in output tables
- `experiment`: (optional) label to use for the "Experiment" column in outputs; useful for downstream data analysis
- `interval_beds`: (optional) array of `.bed` files to use for creating the `jaccard_table` output
- `interval_bed_labels`: (optional) array of labels for the `.bed` files to use in the output tables

### Outputs

- `whatshap_eval`: the evaluation output containing error statistics
- `switch_errors_bed`: a bed file containing locations of switch errors
- `longest_blocks`: a tsv containing the variants in the longest phased block in each contig
- `jaccard_table`: if `.bed` files were provided, a table of the Jaccard index of the `switch_errors_bed` and  its intersection
with the input `.bed` files; this statistic reflects what percentage of the switch errors occured in each of the `.bed` files.

### Tool Options

Some options can be toggled in the `WhatsHapCompare` task for different behaviors. See the 
[tool documentation](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-compare) for more details. Included
options are:

- `only_snvs`: (default = `false`) If true, only use SNVs in phasing comparison