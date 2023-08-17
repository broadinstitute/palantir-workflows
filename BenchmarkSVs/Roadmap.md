# Development Roadmap

This file contains a "to do" list for planned changes and improvements to the WDL and plotting scripts.

## BenchmarkSVs WDL

- [x] Add interval list overlap stats for breakpoints (with padding) of SVs rather than full spanning interval.
- [x] Collect QUAL scores for downstream analysis (as some tools use interpretable measurement).
- [ ] Fix bed file "off by one" error -- currently use bed and VCF coordinates together, so need to shift bed coords by one.
- [x] Add optional bed file for global restriction of variants.
- [ ] Update `query` statements to filter out entries without required INFO fields.

## SVisualizer

- [x] Allow for basic counts plots to filter by interval overlap percents, lengths, etc. (More control).
- [x] Add plots for SVLEN and QUAL distributions.
- [x] Add more interval list overlap sliders (esp. to Basic Wittyer Tab), and for breakpoint overlaps after collecting data.
- [x] Add `ALL` for SVTYPE option in Adv Wittyer Plots.
- [x] Allow user to control fixed vs dynamic axes in Prec/Recall plots.
- [ ] Update HWE plots to use density rather than one dot per variant to prevent crashing with large files.
- [x] Add toggle for Truvari plots to force GT match. (Added GT Concordance plot instead.)

## Backlog

Some tasks that are lower priority at the moment but might get picked up in the future. 

- [ ] QC: Add optional .ped file for collecting Mendelian concordance stats.
- [ ] Add CNV-aware methods for collecting GT info (since CNVs have variable ploidy -- use CN field?).
- [ ] Allow for custom site-level padding around breakpoints when collecting intersection stats.


## CHANGELOG

### v0.8

- Fixed bug in `CleanSVs.wdl` to preserve existing abstract alleles in file. Also convert BND notation using `[` or `]` to abstract.
Also try to infer SVTYPE from allele when UNK (e.g. from GATK-SV). User can now toggle whether to write BND variants or filter them out.
- Record user comment for submission in `gather_terra_data` README file.
- Update docker version to use bcftools 1.18 to avoid an old bug.

### v0.7

- Added column for Terra workflow id when using `gather_terra_data` script to help disambiguate data across runs when using
similar Experiment labels between runs.
- Minor changes to the main WDL.

### v0.6

- Added auto-generated `README.txt` to the directory created by `gather_terra_data.py` script.
- Fixed `CleanSVs.wdl` to remove records where len(REF) - len(ALT) is not above threshold to avoid weird Truvari annotations
if the sequences are large but only have small variants. Also switch negative length to positive (e.g. from sniffles2).
- Fixed bug in `BenchmarkSVs.wdl` with compression when subsetting to an evaluation region.
- Allow user to manually fix order of Experiment categories for coloring plots in SVisualizer.
- Added GT Concordance plots in Truvari benchmarking tab (below main plot), which represents "precision" for genotypes.

### v0.5

- Added a `CleanSVs.wdl` script which can take a VCF with large INDELs (including sequences) and clean to optionally split MA sites, 
add annotations/remove short variants, and replace with abstract alleles/convert missing `.` to `PASS` filters (for Wittyer)
- Small bug fixes/tweaks to the `BenchmarkSVs.wdl` workflow
- Improved logic for rendering HWE plots to allow exclusion of second copy if requested at top of script (for case where comparing
single sample VCFs against a panel)
- Refactored SVisualizer code to plot Counts tab plots as averages with error bars when multiple Base/Comp experiment sample labels exist
- Now collect `QUAL` stats and plot distributions (plus `SVLEN` distributions)
- Added some extra labels for `Experiment` in output to allow for plotting averages across different groups in plots


### v0.4

- Allow for `evaluation_bed` input along with `evaluation_pct` to restrict to variants intersecting a percentage of eval regions globally
- Added toggle to WDL for splitting MA sites to ensure only biallelic sites are used downstream
- Added QC collection for baseline VCF as well, with experiment label separating the two in combined output
- Added overlap stats from `bed_regions` files into QC output files as well (previously inferred in plotting script w/ error-prone logic)
- Changed the Truvari code to fix bug where sample pair stats were not being properly computed
- Made a few tweaks to SVisualizer script (e.g. ensuring categorical orders for binning in counts tab, etc)

### v0.3

- Added a `combined_files.tar.gz` output to capture all top level outputs into one archive containing directory structure
- Updated `gather_terra_data.py` to work at level of submission instead of workflow; combined data across runs inside submission
- Added unified interface for filtering based on interval data, FILTERs, etc including toggle for fixed or dynamic axes on prec/recall
- Added buttons for fixed vs dynamic axes when appropriate
- Allow `COVARIATE_X` input at top of notebook for allowing custom covariate to measure prec/recall against, e.g. coverage (user provided)
- Add extra hover data for scatter plots
- Fix a few bugs / cover some edge cases

### v0.2

- Added breakpoint overlap stats collection to pipeline
- Fixed some title issues for HWE plots
- Implemented `ALL` type logic for Adv Wittyer Plots sidebar plugin

### v0.1

- Initial "release"