# Joint Calling

GATK's joint genotyping pipelines take their inputs as a *sample map*: a two-column, tab-delimited file of
sample name and GVCF path, one sample per line. Assembling that file by hand for a large cohort is tedious, and
in Terra the GVCF paths usually live in a data table column rather than in a file. This directory contains a
small workflow that takes two parallel arrays and writes the sample map for you.

## CreateSampleMapFile

Defined in [CreateSampleMap.wdl](CreateSampleMap.wdl). Note the workflow name is `CreateSampleMapFile`, while
the task inside it is named `CreateSampleMap`.

### Summary

A single task running inline Python. It zips `sample_names` with `gvcfs` and writes one
`<sample_name>\t<gvcf_path>` line per sample to a file named `inputs.list`. If the two arrays are not the same
length, the script prints an error and exits non-zero, failing the workflow.

The GVCF inputs are declared `localization_optional: true` in `parameter_meta`, so the GVCFs are never
downloaded — only their paths are written into the map. This keeps the task cheap regardless of cohort size, and
means the paths in the output map are the original bucket paths.

### Inputs

* `Array[String] sample_names`: sample names, in the same order as `gvcfs`
* `Array[File] gvcfs`: GVCF paths, in the same order as `sample_names`. Not localized (see above), so these
  should be paths readable by whatever downstream tool consumes the map

### Outputs

* `File sample_map_file`: the generated `inputs.list`, a headerless tab-delimited file with one
  `sample_name<TAB>gvcf_path` row per sample

### Notes

* Runtime is fixed in the task and not exposed as workflow inputs: `python:3.6`, 1 CPU, 2 GB memory,
  `local-disk 50 HDD`. There are no preemptible attempts configured.
* The output filename is always `inputs.list`; it is not configurable.
* This workflow is not currently registered in [`.dockstore.yml`](../.dockstore.yml) and has no test entry in
  [`test/watt_config.yml`](../test/watt_config.yml).

## CreateSampleMap_template.json

[`CreateSampleMap_template.json`](CreateSampleMap_template.json) is a minimal, filled-in example of the inputs
JSON, meant to be copied and edited rather than used as-is. It shows the two required keys, fully qualified with
the workflow name:

```json
{
  "CreateSampleMapFile.gvcfs": [
    "gs://bucket/file_HG00096.g.vcf.gz",
    "gs://bucket/file_HG00419.g.vcf.gz"
  ],
  "CreateSampleMapFile.sample_names": [
    "HG00096",
    "HG00419"
  ]
}
```

The two arrays are positional: `sample_names[i]` is paired with `gvcfs[i]`, so the ordering must match. Replace
the `gs://bucket/...` placeholders with your real GVCF paths and the sample names with the corresponding sample
identifiers, keeping the array lengths equal.
