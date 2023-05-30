# Development Roadmap

This file contains a "to do" list for planned changes and improvements to the WDL and plotting scripts.

## BenchmarkSVs WDL

- [ ] Add interval list overlap stats for breakpoints (with padding) of SVs rather than full spanning interval.
- [ ] QC: Add optional .ped file for collecting Mendelian concordance stats.
- [ ] Collect QUAL scores for downstream analysis (as some tools use interpretable measurement).
- [ ] Add CNV-aware methods for collecting GT info (since CNVs have variable ploidy -- use CN field?).
- [ ] Fix bed file "off by one" error -- currently use bed and VCF coordinates together, so need to shift bed coords by one.

## SVisualizer

- [ ] Allow for basic counts plots to filter by interval overlap percents, lengths, etc. (More control).
- [ ] Add plots for SVLEN distributions.
- [ ] Add more interval list overlap sliders (esp. to Basic Wittyer Tab), and for breakpoint overlaps after collecting data.
- [ ] Add `ALL` for SVTYPE option in Adv Wittyer Plots.
- [ ] Allow user to control fixed vs dynamic axes in Prec/Recall plots.