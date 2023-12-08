# Interval Files

This document contains information about popular interval files (.bed) used in restricting analysis to subsets of the
genome, particularly for benchmarking. These are primarily taken fron NIST's GIAB Genome Stratifications latest (v3.1) 
release [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/), 
which contains many more files than included here with in-depth documentation. This page is meant to be a summary of 
"greatest hits" from that site, along with a few others. Note this does not contain the high confidence interval lists
for specific GIAB samples, which are mirrored here: `gs://broad-dsde-methods-hydro-gen-truth-data-public/NIST/GIAB/v4/` .

Some interval files have "interesting overlap" statistics to convey at a glance how much redundancy there may be between
two interval files that have conceptual similarity. The appendix at the end goes into some detail as to how these 
were computed.


## Table of Contents

1. [Interval File Summaries](#interval-file-summaries)
   1. [Reference Specific](#reference-specific)
      1. [HighGC](#highgc)
      2. [LowGC](#lowgc)
      3. [LowMappability](#lowmappability)
      4. [Homopolymers](#homopolymers)
      5. [LowComplexityRegion](#lowcomplexityregion)
      6. [SegDup](#segdup)
      7. [AllDifficult](#alldifficult)
   2. [Function Specific](#function-specific)
      1. [Exome](#exome)
      2. [CodingRegions](#codingregions)
      3. [BadPromoters](#badpromoters)
      4. [CMRG](#cmrg)
   3. [Sequence Process Specific](#sequence-process-specific)
      1. [UGHiConf](#ughiconf)
2. [Using on Terra](#using-on-terra)
3. [Appendix on Statistics](#appendix-on-statistics)


## Interval File Summaries

### Reference Specific

These are interval files which are intrinsic to fixed reference genome sequences (hg38 for this version). Each has three
mirror locations: one bed, one interval_list with hg38 header, and one interval_list with hg38 minus alt contigs for 
its sequence dictionary (no-alt). Full lists of each category are collected in the [Using on Terra](#using-on-terra) section.

#### HighGC
* Summary: Contains regions in reference consisting of 100bp windows with more than 85% GC bases. Note the `slop50` in
the file name references adding 50 bases of padding on either side of the window using 
[bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/slop.html), so the list is the union of 200-length
intervals with interior 100bps having high GC content. Note other thresholds than 85% are available from NIST.
* Bucket Mirrors: 
  * bed: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_gc85_slop50.bed`
  * interval_list: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_gc85_slop50.interval_list`
  * interval_list (no-alt):  `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_gc85_slop50.interval_list`
* Original File: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/GCcontent/GRCh38_gc85_slop50.bed.gz)
* Original Documentation: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/GCcontent/GRCh38-GCcontent-README.md)
* Percent of Genome: 0.01%
* Interesting Overlaps:
  * % of HighGC in LowComplexityRegion: 28% (Jaccard = 0.2%)

#### LowGC
* Summary: Same as HighGC, but with inner windows containing less than 15% GC bases.
* Bucket Mirrors: 
  * bed: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_lt15_slop50.bed`
  * interval_list: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_lt15_slop50.interval_list`
  * interval_list (no-alt): `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_lt15_slop50.interval_list`
* Original File: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/GCcontent/GRCh38_gc15_slop50.bed.gz)
* Original Documentation: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/GCcontent/GRCh38-GCcontent-README.md)
* Percent of Genome: 0.2%
* Interesting Overlaps:
  * % of LowGC in LowComplexityRegion: 60% (Jaccard = 4.5%)

#### LowMappability
* Summary: Contains regions in reference which are homologous to other regions. This is achieved using two sets of 
"stringency" parameters: windows of length >= 100 with multiple alignments having up to 2 mismatch/1 indel < 15 bp 
(low stringency), and regions of length >= 250 with exact alignment matches (high stringency). Note that while 
conceptually one might expect the latter to be a subset of regions in the former, the latter list has a few more sites
than seen in the former (0.03% of high stringency not seen in low). Refinements to either subset is available from NIST. 
* Bucket Mirrors: 
  * bed: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_lowmappabilityall.bed`
  * interval_list: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_lowmappabilityall.interval_list`
  * interval_list (no-alt): `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_lowmappabilityall.interval_list`
* Original File: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/mappability/GRCh38_lowmappabilityall.bed.gz)
* Original Documentation: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/mappability/GRCh38-mappability-README.md)
* Percent of Genome: 8%
* Interesting Overlaps:
    * % of LowMappability in Homopolymers: 2% (Jaccard = 1.5%)
    * % of LowMappability in LowComplexityRegion: 1.5% (Jaccard = 1%)
    * % of LowMappability in SegDup: 44% (Jaccard = 36%)
    * % of LowMappability in AllDifficult: 100% (contained; Jaccard = 39%)

#### Homopolymers
* Summary: Contains regions in reference which are homopolymers (repeated base) of length >= 6 and imperfect homopolymers
  (homopolymers interrupted by 1 different base) of length >= 10, with 5 bases added on both sides for padding.
* Bucket Mirrors: 
  * bed: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed`
  * interval_list: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.interval_list`
  * interval_list (no-alt): `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.interval_list`
* Original File: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/LowComplexity/GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz)
* Original Documentation: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/LowComplexity/GRCh38-LowComplexity-README.md)
* Percent of Genome: 2.7%
* Interesting Overlaps:
  * % of Homopolymers in LowMappability: 6% (Jaccard = 1.5%)
  * % of Homopolymers in LowComplexityRegion: 37% (Jaccard = 26%)
  * % of Homopolymers in SegDup: 5.6% (Jaccard = 2%)
  * % of Homopolymers in AllDifficult: 100% (contained; Jaccard = 13%)

#### LowComplexityRegion
* Summary: Contains regions of "low complexity", using the symmetric DUST algorithm. 
* Bucket Mirrors: 
  * bed: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/LCRFromHengHg38.bed`
  * interval_list: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/LCRFromHengHg38.interval_list`
  * interval_list (no-alt): `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_LCRFromHeng.interval_list`
* Original File: [here](https://github.com/broadinstitute/hydro.gen/blob/main/Data/LCRFromHengHg38.bed) 
* Original Documentation: [here](https://github.com/broadinstitute/hydro.gen/tree/main/Data)
* Percent of Genome: 2%
* Interesting Overlaps:
  * % of LowComplexityRegion in LowMappability: 6% (Jaccard = 1%)
  * % of LowComplexityRegion in Homopolymers: 47% (Jaccard = 26%)
  * % of LowComplexityRegion in AllDifficult: 93.7% (Jaccard = 9.5%)

#### SegDup
* Summary: Contains regions of segmental duplications, i.e. regions where the reference aligns to itself. 
* Bucket Mirrors:  
  * bed: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_segdups.bed`
  * interval_list: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_segdups.interval_list`
  * interval_list (no-alt): `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_segdups.interval_list`
* Original File: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/SegmentalDuplications/GRCh38_segdups.bed.gz)
* Original Documentation: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/SegmentalDuplications/GRCh38-SegDups-README.md)
* Percent of Genome: 5.4%
* Interesting Overlaps:
  * % of SegDup in LowMappability: 66% (Jaccard = 36%)
  * % of SegDup in Homopolymers: 2.8% (Jaccard = 2%)
  * % of SegDup in AllDifficult: 100% (contained; Jaccard = 26%)

#### AllDifficult
* Summary: Contains regions pooled from various "difficult" regions, including all tandem repeats, homopolymers of length >= 6,
segmental duplications, and more. Note this includes a few "function specific" interval lists like BadPromoters.
* Bucket Mirrors: 
  * bed: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_alldifficultregions.bed`
  * interval_list: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_alldifficultregions.interval_list`
  * interval_list (no-alt): `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_alldifficultregions.interval_list`
* Original File: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/Union/GRCh38_alldifficultregions.bed.gz)
* Original Documentation: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/Union/GRCh38-Union-README.md)
* Percent of Genome: 21%
* Interesting Overlaps:
  * % of AllDifficult in LowMappability: 39% (Jaccard = 39%)
  * % of AllDifficult in Homopolymers: 13% (Jaccard = 13%)
  * % of AllDifficult in LowComplexityRegion: 9.6% (Jaccard = 9.5%)
  * % of AllDifficult in SegDup: 26% (Jaccard = 26%)
  * % of AllDifficult in UGHiConf: 74% (Jaccard = 16%)


### Function Specific

These are interval files which are related to inferred biological functioning of the underlying sequences.

#### Exome
* Summary: Contains exome with additional padding (50bp) and a few other modifications. 
* Bucket Mirrors: 
  * bed: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/exome_evaluation_regions.v1.bed`
  * interval_list: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/exome_evaluation_regions.v1.interval_list`
  * interval_list (no-alt): `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/exome_evaluation_regions.v1-noalt.interval_list`
* Original File: Same as bucket.
* Original Documentation: Script to generate found [here](https://github.com/broadinstitute/warp/blob/ldgauthier-patch-splitintervals/scripts/generate_hg38_exome_calling_intervals.v1.1.sh)
* Percent of Genome: 1.14%
* Interesting Overlaps:
  * % of Exome in AllDifficult: 28% (Jaccard = 1.5%)
  * % of Exome in CodingRegions: 95% (Jaccard = 94%)
  * % of Exome in UGHiConf: 97% (Jaccard = 1.2%)
  * % of Exome in CMRG: 1.8%

#### CodingRegions
* Summary: Contains regions strictly formed from coding regions, as identified by RefSeq. In essence, this is conceptually 
a stricter exome interval file, though there is not a strict containment.
* Bucket Mirrors: 
  * bed: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_refseq_cds.bed`
  * interval_list: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_refseq_cds.interval_list`
  * interval_list (no-alt): `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_refseq_cds.interval_list`
* Original File: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/FunctionalRegions/GRCh38_refseq_cds.bed.gz)
* Original Documentation: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/FunctionalRegions/GRCh38-FunctionalRegions-README.md)
* Percent of Genome: 1.10%
* Interesting Overlaps:
  * % of CodingRegions in Exome: 98.6% (Jaccard = 94%)

#### BadPromoters
* Summary: Contains regions where transcription site or first exons have systematically low coverage.
* Bucket Mirrors: 
  * bed: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_BadPromoters.bed`
  * interval_list: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_BadPromoters.interval_list`
  * interval_list (no-alt): `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_BadPromoters.interval_list` 
* Original File: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/FunctionalTechnicallyDifficultRegions/GRCh38_BadPromoters.bed.gz)
* Original Documentation: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/FunctionalTechnicallyDifficultRegions/GRCh38-FunctionalTechnicallyDifficult-README.md)
* Percent of Genome: 0.006%
* Interesting Overlaps:
  * % of BadPromoters in HighGC: 25% (Jaccard = 9%)
  * % of BadPromoters in LowGC: 0% (Jaccard = 0%)
  * % of BadPromoters in LowMappability: 2% (Jaccard = 0.002%)
  * % of BadPromoters in Homopolymers: 7% (Jaccard = 0.015%)
  * % of BadPromoters in LowComplexity: 27% (Jaccard = 0.08%)
  * % of BadPromoters in SegDup: 4% (Jaccard = 0.0055%)
  * % of BadPromoters in Exome: 6.4% (Jaccard = 0.035%)
  * % of BadPromoters in CodingRegions: 6.5% (Jaccard = 0.037%)

#### CMRG
* Summary: A list of Challenging Medically Relevant Genes sequenced specifically for HG002. Note this also pairs with a separate truth VCF for HG002.
* Bucket Mirrors:
  * bed: `gs://broad-dsde-methods-hydro-gen-truth-data-public/NIST/GIAB_CMRG/v1.00/HG002_GRCh38_CMRG_smallvar_v1.00.bed`
  * interval_list: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/HG002_GRCh38_CMRG_smallvar_v1.00.interval_list`
  * interval_list (no-alt): `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/HG002_GRCh38-noalt_CMRG_smallvar_v1.00.interval_list`
  * truth vcf: `gs://broad-dsde-methods-hydro-gen-truth-data-public/NIST/GIAB_CMRG/v1.00/HG002_GRCh38_CMRG_smallvar_v1.00.broad-header.vcf.gz`
* Original File: [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.bed)
* Original Documentation: [here](https://data.nist.gov/od/id/mds2-2475)
* Percent of Genome: 0.38%
* Interesting Overlaps:
  * % of CMRG in Exome: 5.4%


### Sequence Process Specific

These are interval lists that are relevant only for analyzing data from specific sequencing methods.

#### UGHiConf
* Summary: Contains regions defined by Ultima Genomics to be of "high confidence", i.e. removing the union of various
"difficult" regions for Ultima sequencers. The full list can be found in the original documentation, but this includes
homopolymers of length >= 11 (+4 bases for padding), GC < 5%, and 50bp windows with highly variable coverage between
samples.
* Bucket Mirrors: 
  * bed: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/ug_hcr.bed`
  * interval_list: `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/ug_hcr.interval_list`
  * interval_list (no-alt): `gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/ug_hcr-noalt.interval_list`
* Original File: Same as bucket.
* Original Documentation: [here](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/ultima_genomics/interval_lists/README.txt)
* Percent of Genome: 89.5%
* Interesting Overlaps:
  * % of UGHiConf in AllDifficult: 17% (Jaccard = 16%)

## Using on Terra

Here is a formatted list of the above interval files along with a list of names in the same order, so you can easily
use any subset of them when running a benchmarking tool on Terra. Depending on your workflow, you might require bed
or interval list files, so we provide collected lists of both.

```
strat_intervals_bed = [
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_gc85_slop50.bed",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_lt15_slop50.bed",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_lowmappabilityall.bed",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/LCRFromHengHg38.bed",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_segdups.bed",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_alldifficultregions.bed",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/exome_evaluation_regions.v1.bed",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_refseq_cds.bed",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_BadPromoters.bed",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/NIST/GIAB_CMRG/v1.00/HG002_GRCh38_CMRG_smallvar_v1.00.bed",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/ug_hcr.bed"
]

# Interval Lists w/ full hg38 sequence dictionary
strat_interval_lists = [
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_gc85_slop50.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_lt15_slop50.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_lowmappabilityall.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/LCRFromHengHg38.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_segdups.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_alldifficultregions.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/exome_evaluation_regions.v1.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_refseq_cds.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38_BadPromoters.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/HG002_GRCh38_CMRG_smallvar_v1.00.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/ug_hcr.interval_list"
]

# Interval Lists w/ hg38 sequence dictionary without alt contigs
strat_interval_lists_no_alt = [
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_gc85_slop50.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_lt15_slop50.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_lowmappabilityall.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_LCRFromHeng.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_segdups.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_alldifficultregions.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/exome_evaluation_regions.v1-noalt.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_refseq_cds.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/GRCh38-noalt_BadPromoters.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/HG002_GRCh38-noalt_CMRG_smallvar_v1.00.interval_list",
"gs://broad-dsde-methods-hydro-gen-truth-data-public/IntervalFiles/ug_hcr-noalt.interval_list"
]

strat_labels = [
"HighGC",
"LowGC",
"LowMappability",
"Homopolymers",
"LowComplexityRegion",
"SegDup",
"AllDifficult",
"Exome",
"CodingRegions",
"BadPromoters",
"CMRG",
"UGHiConf",
]
```

## Appendix on Statistics

Above, we compute the amount of redundancy between two interval lists using the [Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index).
Recall for two sets A and B, this is given by the size of A intersect B divided by the size of A union B. Note this is bounded
by 0 and 1. It is equal to 1 if and only if the two sets are identical, and equal to 0 if and only if they are disjoint.

We compute this statistic on a few "interesting" pairs and report it above, along with asymmetric redundancy calculations
measuring what percentage of A lies in B. This is the special case of the Jaccard index of A intersect B with A. Note also
the percent of genome covered can be computed using a Jaccard index of your bed file A with a "full" bed file covering
all contigs of your genome. A script to make a "full" bed from a reference dict is included in the `create_full_bed.sh`
file. In the above calculations, a "full" bed was created from hg38 for just the autosomes and X/Y to compute what 
percent of the genome a given interval file spanned. Note that this includes bases labeled as `N` in the reference.
To get a percent of the non-`N` bases, you can subset hg38 to exclude them (and then run the same `overlap_stats.sh`
script described below).

The script to compute all three of these (the usual Jaccard index, percent of A in B, and percent of B in A) for two 
input bed files A and B is provided in the `overlap_stats.sh` file for reproducibility of the above code and extension 
beyond the provided cases. You must have the `hg38-noalt.genome` or similar file for your reference in the same directory 
when running this script for it to function properly, which is an ordered list of contigs with sequence length. Overlap
stats are only computed over the "no-alt" version of hg38.
