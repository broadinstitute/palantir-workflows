# Fingerprints

Picard's fingerprinting tools genotype a sample at a small, fixed panel of sites (a "haplotype map") so that
files can later be checked for sample identity. Before relying on those fingerprints, it is useful to know how
informative and how self-consistent a given fingerprint actually is. This directory contains a thin wrapper that
runs Picard's `CalculateFingerprintMetrics` over a fingerprint VCF and returns the resulting metrics file.

## GetFingerprintMetrics

### Summary

A single-task workflow. It runs

```
java -jar <picard_jar> CalculateFingerprintMetrics INPUT=<fingerprint_vcf> OUTPUT=<output_prefix>.fingerprint_metrics HAPLOTYPE_MAP=<haplotype_database>
```

and returns the metrics file produced by Picard. There is no pre- or post-processing done by the WDL.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Fingerprints/GetFingerprintMetrics.html) · [open locally](../docs/viz/Fingerprints/GetFingerprintMetrics.html)

### Inputs

* `File fingerprint_vcf`: the fingerprint VCF to compute metrics for; passed to Picard as `INPUT`
* `File picard_jar`: the Picard jar to run. This is a required input — the task runs on the stock `openjdk:8`
  image, which contains no Picard installation, so the jar must be supplied and must be a build that includes
  the `CalculateFingerprintMetrics` tool
* `File haplotype_database`: the haplotype map defining the fingerprinting sites; passed to Picard as `HAPLOTYPE_MAP`.
  This must be the same haplotype map the fingerprint VCF was generated against
* `String output_prefix`: prefix for the output file, which is named `<output_prefix>.fingerprint_metrics`

### Outputs

* `File fingerprint_metrics`: the `<output_prefix>.fingerprint_metrics` file written by Picard

### Notes

* Runtime is fixed in the task and not exposed as workflow inputs: `openjdk:8`, 1 CPU, 2 GB memory,
  `local-disk 100 LOCAL`, 3 preemptible attempts.
* Haplotype maps suitable for the `haplotype_database` input can be produced with
  [BuildHaplotypeMap](../HaplotypeMap/README.md).
* To compare fingerprints *between* files rather than characterize a single one, see
  [MatchFingerprints](../Utilities/WDLs/README.md#matchfingerprints), which wraps Picard's
  `CrosscheckFingerprints` and takes the same kind of haplotype map file.
* This workflow is not currently registered in [`.dockstore.yml`](../.dockstore.yml) and has no test entry in
  [`test/watt_config.yml`](../test/watt_config.yml).
