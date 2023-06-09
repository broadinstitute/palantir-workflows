# Development Roadmap

This file contains a "to do" list for planned changes and improvements to the WDL and plotting scripts.

## BenchmarkSVs WDL

- [x] Add interval list overlap stats for breakpoints (with padding) of SVs rather than full spanning interval.
- [ ] QC: Add optional .ped file for collecting Mendelian concordance stats.
- [ ] Collect QUAL scores for downstream analysis (as some tools use interpretable measurement).
- [ ] Add CNV-aware methods for collecting GT info (since CNVs have variable ploidy -- use CN field?).
- [ ] Fix bed file "off by one" error -- currently use bed and VCF coordinates together, so need to shift bed coords by one.

## SVisualizer

- [x] Allow for basic counts plots to filter by interval overlap percents, lengths, etc. (More control).
- [ ] Add plots for SVLEN distributions.
- [x] Add more interval list overlap sliders (esp. to Basic Wittyer Tab), and for breakpoint overlaps after collecting data.
- [x] Add `ALL` for SVTYPE option in Adv Wittyer Plots.
- [x] Allow user to control fixed vs dynamic axes in Prec/Recall plots.
- [ ] Update HWE plots to use density rather than one dot per variant.
- [ ] Add toggle for Truvari plots to force GT match.

## CHANGELOG

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