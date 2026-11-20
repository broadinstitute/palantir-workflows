# ScoreBGE

Polygenic scoring for BGE (Blended Genome Exome) data, where each sample has both a whole exome sequencing (WES) GVCF
and an imputed whole genome (WGS) VCF. This directory contains both the WDL workflow and the python package it runs.

See the [PRS README](../README.md) for how this fits into the wider PRS pipeline. `ScoreBGE` is called by
[ScoringPart.wdl](../ScoringPart.wdl) (when `use_bge_scoring` is true) and by [PCARE.wdl](../PCARE.wdl), and is
registered on Dockstore as **ScoreBGE** (see [../../.dockstore.yml](../../.dockstore.yml)).

## Scoring logic

The scorer walks the weights file in reference-dictionary order and, for each weight:

1. **WES GVCF first.** It looks for a record (variant record or reference block) covering the site. A sample is scored
   from the GVCF only if its `GQ` at that site is at least the GQ threshold (default 30) and its genotype is fully
   called; no-calls and half-calls count as low quality and are skipped. The site score is the number of copies of the
   effect allele times the weight. With `score_haploid_as_diploid`, a single-allele genotype (e.g. on chrX) is scored
   as if it were homozygous.
2. **Imputed WGS VCF second, only where the exome did not score.** For each sample, sites already scored from the GVCF
   are skipped. Remaining sites are scored from the `DS` dosage field (flipped to `2 - DS` when the effect allele is
   the reference allele); the task errors if `DS` is absent. Only biallelic records are supported.
3. The final score for a sample is the sum of its GVCF and VCF scores.

Sample names must be identical (and in the same order) between the GVCF and the VCF, or an explicit `--sample-names`
list must be given to both passes. The scorer logs per-source metrics (min/max sites scored, low quality sites, sites
not found in the VCF) to stdout.

## ScoreBGE (WDL)

Defined in [ScoreBGE.wdl](ScoreBGE.wdl). The workflow is a thin wrapper around the single `ScoreGvcfAndVcf` task,
which runs `python3 /ScoreBGE.py` inside the ScoreBGE docker image.

### Summary

Scores one set of samples against a weights file using both a WES GVCF and an imputed WGS VCF, preferring the exome
genotype wherever it is high quality, and emits per-source and combined plink-style score tables plus lists of the
sites that were scored.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/PRS/ScoreBGE/ScoreBGE.html) · [open locally](../../docs/viz/PRS/ScoreBGE/ScoreBGE.html)

### Inputs

- **File exome_gvcf**: WES GVCF to score.
- **File exome_gvcf_index**: index for `exome_gvcf`.
- **File imputed_wgs_vcf**: imputed WGS VCF to score; must carry `DS` dosages.
- **File imputed_wgs_vcf_index**: index for `imputed_wgs_vcf`.
- **String basename**: prefix for all output file names.
- **File weights**: weights file. The native format is a TSV with header columns `contig`, `position`, `ref`, `alt`,
  `effect_allele`, `weight`.
- **Array[String]? sample_names**: restrict scoring to these samples; if absent, all samples in the files are scored.
- **Boolean score_haploid_as_diploid**: always score haploid genotypes (such as on chrX) as if diploid.
- **Boolean use_emerge_weight_format = false**: read the weights file in the eMERGE format instead of the native
  format. Passed to the script as `--use-emerge-weight-format`.
- **String? score_bge_docker**: override the docker image. Defaults to
  `us.gcr.io/broad-dsde-methods/palantir-workflows-score-bge:palantir-workflows_0480e5e`.
- **File ref_dict**: reference sequence dictionary; the `@SQ` lines define contig ordering for the weights.
- **Int preemptible = 1**: preemptible attempts.

The `ScoreGvcfAndVcf` task additionally exposes `Int? disk_gb` (default: sizes of the GVCF, VCF and weights plus
50 GiB), `Int mem_gb = 4` and `Int cpu = 4`, which are not surfaced as workflow-level inputs.

### Outputs

- **File exome_gvcf_score**: `<basename>.exome_gvcf.score` — `#IID` / `SCORE1_SUM` table of the exome-derived score.
- **File imputed_wgs_vcf_score**: `<basename>.imputed_wgs_vcf.score` — the same for the imputed-VCF-derived score.
- **File score**: `<basename>.score` — the combined score (exome + imputed).
- **File exome_gvcf_sites_scored**: `<basename>.exome_gvcf.sites_scored` — TSV of `site` and the comma-separated list
  of samples scored at that site from the GVCF.
- **File imputed_wgs_vcf_sites_scored**: `<basename>.imputed_wgs_vcf.sites_scored` — the same for the imputed VCF.
- **File any_source_any_sample_sites_scored**: `<basename>.any_source_any_sample.sites_scored` — a plain, plink
  compatible list of `contig:position:ref:alt` IDs for every site scored from either source in any sample. This is
  what `ScoringPart.wdl` uses as the "sites scored" list when comparing against the sites used in model training.

## ScoreBGE.py (python package)

[ScoreBGE.py](ScoreBGE.py) provides the `BGEScorer` class and a command-line entry point. Dependencies are pinned in
[requirements.txt](requirements.txt): `pysam==0.20.0`, `pandas==1.3.4`, `numpy==1.21.4`.

### Command line

```
python3 ScoreBGE.py \
  --weights <weights.tsv> \
  --gvcf <exome.g.vcf.gz> \
  --vcf <imputed_wgs.vcf.gz> \
  --ref-dict <reference.dict> \
  --basename <output_basename> \
  [--sample-names SAMPLE [SAMPLE ...]] \
  [--score-haploid-as-diploid]
```

The `--weights`, `--gvcf`, `--vcf`, `--ref-dict` and `--basename` arguments are required. The CLI runs
`score_wes_gvcf`, then `score_wgs_vcf`, then `write_output`.

**Note:** the WDL also passes `--use-emerge-weight-format` when `use_emerge_weight_format` is set, but that flag is
not present in the version of `ScoreBGE.py` committed here. The pinned docker image
(`palantir-workflows-score-bge:palantir-workflows_0480e5e`) is built from a version of the script that supports it. If
you rebuild the image from this directory, check that `use_emerge_weight_format` is still supported before using it.

### API

- `BGEScorer(ref_dict_path, prs_weights_path, output_basename, score_haploid_as_diploid)`: reads the reference
  dictionary and the weights (sorting weights by contig order, then position).
- `score_wes_gvcf(gvcf_path, sample_names=None, site_gq_threshold=30)`: scores the exome GVCF and writes
  `<basename>.exome_gvcf.sites_scored`. Must be called before `score_wgs_vcf` unless `allow_wgs_vcf_only` is used.
  Raises if either pass has already been run on this object.
- `score_wgs_vcf(wgs_vcf_path, sample_names=None, allow_wgs_vcf_only=False)`: scores the imputed WGS VCF at sites not
  already scored per-sample from the GVCF, and writes `<basename>.imputed_wgs_vcf.sites_scored`. Raises if sample
  names do not match those from the GVCF pass.
- `write_output(allow_single_source_scoring=False)`: writes the two per-source `.score` tables, the combined
  `.score` table and `<basename>.any_source_any_sample.sites_scored`. Raises if either source was not scored unless
  `allow_single_source_scoring` is set.

### Docker

[Dockerfile](Dockerfile) is `FROM us.gcr.io/broad-dsde-methods/python-data-slim-pysam:v1.0` and adds `ScoreBGE.py` at
`/ScoreBGE.py`. Build and push it with the shared script:

```
../build_push_docker.sh --directory ScoreBGE --image-version-tag <tag>
```

Note that this tags the image as `us.gcr.io/broad-dsde-methods/ScoreBGE:<tag>`, whereas the image referenced by the
WDL is named `palantir-workflows-score-bge`; set `score_bge_docker` explicitly if you build your own.

### Tests

pytest tests live in [tests/test_ScoreBGE.py](tests/test_ScoreBGE.py) with fixtures under
[tests/resources](tests/resources) covering WES GVCF variant records, WES GVCF reference blocks, WGS VCF only, and the
combined case, each with expected score and sites-scored files.
