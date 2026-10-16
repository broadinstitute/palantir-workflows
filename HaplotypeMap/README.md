# HaplotypeMap

Picard's fingerprinting tools need a *haplotype map*: a file listing a panel of sites, with a sequence dictionary
header, that is used to genotype and compare samples. A good panel is made of common sites that are close to
independent of each other, so this directory contains a workflow that derives such a panel from a genotyped
multi-sample VCF by filtering to well-called common SNPs and LD-pruning them. It also contains the `Dockerfile`
for the small Python image used by one of the steps.

## BuildHapMap

Defined in [BuildHaplotypeMap.wdl](BuildHaplotypeMap.wdl). Note the workflow name is `BuildHapMap`, not
`BuildHaplotypeMap`, so inputs must be prefixed `BuildHapMap.` even though it is registered on Dockstore under
the name `BuildHaplotypeMap` (see [`.dockstore.yml`](../.dockstore.yml)).

### Summary

Five tasks run in sequence:

1. **`vcftools`** — drops indels and any variant not called in every sample
   (`--remove-indels --max-missing-count <max_missing>`), writing a recoded VCF.
2. **`bcftools`** — `bcftools annotate --set-id '%CHROM:%POS'`, so every site has an ID of the form `chr:pos`.
3. **`plink`** — `plink1.9 --snps-only --biallelic-only --maf <min_maf> --indep-pairwise <prune_window> <prune_slide> <prune_cutoff>`,
   producing a `.prune.in` list of sites that survive LD and minor allele frequency pruning.
4. **`reformat`** — inline Python that writes the haplotype map: the provided sequence dictionary is copied in as
   the header, followed by a `#CHROM POS ID REF ALT INFO` column header and one row per pruned site. Sites are
   further filtered here: multi-base and `*` alternate alleles are skipped. The `INFO` column is filled with the
   first value of the site's `AF` INFO field, so the input VCF must carry `AF`.
5. **`make_vcf`** — Picard `ConvertHaplotypeDatabaseToVcf` turns the haplotype map into a VCF representation.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/HaplotypeMap/BuildHaplotypeMap.html) · [open locally](../docs/viz/HaplotypeMap/BuildHaplotypeMap.html)

### Inputs

* `File input_vcf`: genotyped multi-sample VCF (or GVCF) to derive the panel from. Uncompressed — `vcftools` is
  invoked with `--vcf`. Must have an `AF` INFO field
* `File sequence_dict`: sequence dictionary used verbatim as the haplotype map header; it should match the header
  of the files you intend to fingerprint
* `File reference`: reference FASTA, passed to `ConvertHaplotypeDatabaseToVcf` as `R`
* `File reference_index`: index for `reference`; localized so Picard can find it
* `File reference_dict`: sequence dictionary for `reference`; localized so Picard can find it
* `File picard_jar`: Picard jar used for `ConvertHaplotypeDatabaseToVcf`
* `String output_prefix`: prefix for the two outputs, `<output_prefix>.hapmap.txt` and `<output_prefix>.hapmap.vcf`
* `String? intermediates_prefix`: (task default: `"int"`) prefix for the intermediate files passed between tasks
* `Int? max_missing`: (task default: `0`) maximum number of missing genotypes a site may have; the default keeps
  only sites called in every sample
* `Int? prune_window`: (task default: `50`) `--indep-pairwise` window size
* `Int? prune_slide`: (task default: `5`) `--indep-pairwise` step size
* `Float? prune_cutoff`: (task default: `0.5`) `--indep-pairwise` r^2 threshold
* `Float? min_maf`: (task default: `0.4`) minimum minor allele frequency

Per the header comment in the WDL, higher `prune_window` / `prune_slide` result in more pruning, while lower
`prune_cutoff` / `min_maf` result in less. See the [PLINK LD docs](https://www.cog-genomics.org/plink/2.0/ld)
for details on the pruning parameters.

### Outputs

* `File map`: `<output_prefix>.hapmap.txt`, the haplotype map in Picard's haplotype database format
* `File map_vcf`: `<output_prefix>.hapmap.vcf`, the same panel as a VCF

### Notes

* The `map` output is what you feed to the `haplotype_map` / `haplotype_database` inputs of the fingerprinting
  workflows: [MatchFingerprints](../Utilities/WDLs/README.md#matchfingerprints) and
  [GetFingerprintMetrics](../Fingerprints/README.md).
* Task dockers are pinned per task: `biocontainers/vcftools:v0.1.16-1-deb_cv1`,
  `biocontainers/bcftools:v1.9-1-deb_cv1`, `biocontainers/plink1.9:v1.90b3.45-170113-1-deb_cv1`,
  `docker.io/mollysacks/python:hapmap_builder` (see below), and `broadinstitute/picard` (unpinned tag).
* Runtime is hardcoded in each task and not exposed as workflow inputs; disk requests are large
  (1000–1500 GB HDD) and `make_vcf` asks for 32 GB memory.
* There is no test entry for this workflow in [`test/watt_config.yml`](../test/watt_config.yml).

## Dockerfile

[`Dockerfile`](Dockerfile) builds a small Python image. Nothing in the repo tags it explicitly, but its contents
match what the `reformat` task needs, so it is almost certainly the source for that task's
`docker.io/mollysacks/python:hapmap_builder` image. It is based on `python:3.8.3-slim-buster`, copies the contents of
this directory into `/src`, and `pip install`s `argparse`, `PyVCF`, and `pandas` — the `pandas` and `vcf`
(PyVCF) modules that the inline `reformat` script imports. The script itself lives in the WDL, not in the image,
so nothing under `/src` is actually executed at runtime.

To rebuild and publish:

```bash
cd HaplotypeMap
docker build -t <your-registry>/<your-repo>:<tag> .
docker push <your-registry>/<your-repo>:<tag>
```

Pushing to `mollysacks/python` requires access to that Docker Hub account; in practice you should build under a
registry you control and update the `docker` line in the `reformat` task of
[BuildHaplotypeMap.wdl](BuildHaplotypeMap.wdl) to point at the new image.
