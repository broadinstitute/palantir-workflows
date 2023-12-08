# GLIMPSE Imputation Pipeline

This pipeline performs imputation using the [GLIMPSE2 software](https://github.com/odelaneau/GLIMPSE) by Simone Rubinacci et al. [\[1\]](https://www.biorxiv.org/content/10.1101/2022.11.28.518213v1) [\[2\]](https://www.nature.com/articles/s41588-020-00756-0).

This directory primarily contains the pipeline for GLIMPSE2. It also contains a WDL to run GLIMPSE1 which functions correctly but has not been extensively tested or optimized.

This directory contains the following WDLs:
- [Glimpse2SplitReference](#Glimpse2SplitReference): Split reference panel into chunks and generate binary files to be used in imputation
- [Glimpse2Imputation](#Glimpse2Imputation): Run GLIMPSE2 imputation
- [Glimpse2MergeBatches](#Glimpse2MergeBatches): Combine multiple batches of imputed multi-sample VCFs and recompute AF and INFO annotations
- [Glimpse2ImputationInBatches](#Glimpse2ImputationInBatches): Split inputs into batches, run GLIMPSE2 imputation, and merge the output of the batches back together
- *ReduceAndMergeForGlimpse*: Remove all annotations except GT and PL from single-sample VCFs and merge them into a multi-sample VCF
- *Glimpse1Imputation*: Run GLIMPSE1 imputation (not well-tested or optimized)

## Glimpse2SplitReference

This workflow splits the provided reference panel into a series of chunks that each span a specified region of the genome, and converts the reference panel into a binary representation that can be ingested by the `Glimpse2Imputation` workflow. This significantly speeds up the imputation process since this step only has to be performed once per reference panel instead of once per sample.

**Note on chunk sizes:** The choice of the chunk size is a trade-off between runtime per chunk and number of chunks (and therefore number of parallel jobs to be started). For imputation of single samples we would want to choose a larger minimum window size since each chunk will be processed in a relatively short amount of time and the overhead of starting many jobs for smaller chunks will be considerable. For large batches, however, a smaller minimum window size might be more desireable in order to parallelize more the longer running jobs.

### Input
**Note on paths to the reference panel and genetic maps:** The paths to the reference panel and the genetic maps are constructed by appending a path prefix, the contig name, and a path suffix. For an example path to the reference panel (and index) for chr1 `gs://bucket/to/panel/reference_panel_chr1_merged.bcf(.csi)` the following values would have to be set:
```
reference_panel_prefix = "gs://bucket/to/panel/reference_panel_"
reference_panel_suffix = "_merged.bcf"
reference_panel_index_suffix = ".csi"
contig_name_in_reference_panel = "chr1"
```

- **Array[String] contig_regions**: Defines dinstinct intervals that should be chunked up. Usually, these intervals would consist of whole chromosomes (e.g. `["chr1", "chr2", ...]`). However, for regions like the pseudo-autosomal regions (PAR) of chrX, it might be desireable to split a chromosome into multiple intervals, such as `["chr1", "chr2", ..., "chr22", "chrX:1-2781479", "chrX:2781480-155701382", "chrX:155701383-156030895"]`.
- **Array[String] contig_names_in_reference_panel**: In order to find the correct paths for the reference panel it is necessary to define the contig names that will be inserted between the reference panel perfix and suffix. For the above example this array might look like this: `["chr1", "chr2", ..., "chr22", "chrX", "chrX", "chrX"]`.
- **Array[String] contig_names_in_genetic_maps**: Same as `contig_names_in_reference_panel` but for the genetic maps paths. Following the naming of the genetic maps included in GLIMPSE2 this array would look like this: `["chr1", "chr2", ..., "chr22", "chrX_par1", "chrX", "chrX_par2"]`.
- **String reference_panel_prefix**: See note on paths to the reference panel and genetic maps.
- **String reference_panel_suffix**: See note on paths to the reference panel and genetic maps.
- **String reference_panel_index_suffix**: See note on paths to the reference panel and genetic maps.
- **String genetic_map_path_prefix**: See note on paths to the reference panel and genetic maps.
- **String genetic_map_path_suffix**: See note on paths to the reference panel and genetic maps.
- **Int? seed**: Optional integer seed for the generation of chunks
- **Int? min_window_cm**: Optional minimum window size in [Centimorgan](https://en.wikipedia.org/wiki/Centimorgan). See note on chunk sizes above.
- **Boolean uniform_number_variants = false**: When set to true, each chunk will have approximately the same number of sites while. Each chunk will cover a different (Centimorgan) genetic linkage region with the smallest chunk still being larger than `min_window_cm`.
- **Int preemtible = 1**: Number of preemptible attempts
- **File? monitoring_script**: Optional path to monitoring script. If ommitted, no monitoring will occur and the `split_reference_monitoring** output will not be available.

### Output
- **Array[File] chunks**: One `chunks_contigindex_{CONTIGINDEX}.txt` file per region as defined in the `contig_regions` input that each include the definitions of the chunks in that region.
- **Array[File] split_reference_chunks**: All binary representations of the reference panel for each chunk in each contig region. Each file is named `reference_panel_contigindex_${CONTIGINDEX}_chunkindex_${CHUNKINDEX}` with `CONTIGINDEX` and `CHUNKINDEX` being 4-digit zero-based indices with leading zeros in order to simplify the correct ordering of intervals.
- **Array[String] num_sites**: The number of sites in each chunk.
- **Array[String] num_sites_uniform**: The number of sites in each chunk when using `uniform_number_variants`, otherwise empty.
- **File? monitoring**: A monitoring log if the `monitoring_script` input is set, otherwise null.

## Glimpse2Imputation

This workflow performs the actual imputation using the reference panel chunks generated by GlimpseSplitReference. Each chunk is processed in parallel, subsequently all chunks are being ligated into a single VCF across all regions defined in the `contig_regions` input in GlimpseSplitReference.

The input to this workflow can bei either single-sample or multi-sample VCFs with existing GT and PL annotations or CRAMs, in which case GLIMPSE2 will calculate the PLs using pileup calling.

**Note**: _GLIMPSE2 does not support the input of multiple CRAM files with the same basename when streaming (e.g. `gs://a/file.cram`, `gs://b/file.cram`), due to the way that htslib is implemented. This workflow will check for a potential filename collision and will fail with an error message if such a collision occurs._

### Input
- **File reference_chunks**: File with paths to the `split_reference_chunks` output generated by the GlimpseSplitReference workflow, one file per line. It is necessary to pass them as a list of files because Terra does not like input arrays with 500+ items.
- **File? input_vcf**: VCF to be imputed, single sample or multi sample. Sites not present in the reference panel will be ignored, as are sites present in the reference panel and not in the input VCF. You may only set either this argument or `crams`.
- **File? input_vcf_index**: Index to `input_vcf`.
- **Array[File]? crams**: Array of CRAM files as an alternative to `input_vcf`.
- **Array[File]? cram_indices**: Indices to `crams`. You may only set either this argument or `input_vcf`.
- **Array[String] sample_ids**: Sample IDs to be assigned when using CRAM input.
- **File ref_dict**: By default, GLIMPSE2 only adds the contig names defined by the `contig_regions` input in GlimpseSplitReference to the VCF header's sequence dictionary. In order to make the file compatible with downstream tools the sequence dictionary will be replaced with this input.
- **File? fasta**: Highly recommended when using CRAM input, otherwise the task will attempt to download the reference from a remote server, which can be incredibly slow.
- **File? fasta_index**: Index to `fasta`.
- **String output_basename**: Basename for imputed output VCF.
- **Boolean collect_qc_metrics = true**: If true, runs Hail's [`sample_qc`](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc) method on the imputed VCF and provides the metrics as a TSV table in the `qc_metrics` output.
- **Boolean impute_reference_only_variants = false**: See [GLIMPSE2 documentation](https://odelaneau.github.io/GLIMPSE/docs/documentation/phase/#input-parameters).
- **Boolean call_indels = false**: See [GLIMPSE2 documentation](https://odelaneau.github.io/GLIMPSE/docs/documentation/phase/#bamcram-options-and-filters).
- **Int? n_burnin**: See [GLIMPSE2 documentation](https://odelaneau.github.io/GLIMPSE/docs/documentation/phase/#model-parameters).
- **Int? n_main**: See [GLIMPSE2 documentation](https://odelaneau.github.io/GLIMPSE/docs/documentation/phase/#model-parameters).
- **Int? effective_population_size**: See [GLIMPSE2 documentation](https://odelaneau.github.io/GLIMPSE/docs/documentation/phase/#model-parameters).
- **Int preemtible = 9**: Number of preemptible attempts
- **File? monitoring_script**: Optional path to monitoring script. If ommitted, no monitoring will occur and the `split_reference_monitoring** output will not be available.

### Output

- **File imputed_vcf**: Single imputed VCF that covers all regions defined in the `contig_regions` input in GlimpseSplitReference. The name of the file is the basename of `input_vcf` with `.imputed.vcf.gz` added.
- **File imputed_vcf_index**: Index to `imputed_vcf`.
- **File? qc_metrics**: Output of Hail's [`sample_qc`](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc) method as a TSV table if `collect_qc_metrics` is set, otherwise null.
- **Array[File?] glimpse_phase_monitoring**: A monitoring log for each parallelized chunk if the `monitoring_script` input is set, otherwise null.
- **File? glimpse_ligate_monitoring**: A monitoring log for the ligate task if the `monitoring_script` input is set, otherwise null.

## Glimpse2MergeBatches

This workflow merges multiple batches of imputed multi-sample VCFs into one and recalculates the AF and INFO score annotations based on [this calculation](#glimpse2mergebatches-af-and-info-score-recalculation). If QC metrics for those batches are provided this workflow will also concatenate those metrics into one QC metrics table.

### Input
- **Array[File] imputed_vcf**: GLIMPSE2 output imputed VCFs.
- **Array[File] imputed_vcf_indices**: Indices to `imputed_vcf`.
- **Array[File?]? qc_metrics**: Optional array of qc_metrics of the individual batches to be merged into one `merged_qc_metrics` TSV file. Can otherwise be _null_ or can only contain _null_ values.
- **String output_basename**: Basename for merged output VCF.

### Output
- **File merged_imputed_vcf**: Merged VCF with recalculated AF and INFO annotations.
- **File merged_imputed_vcf_index**: Index to `imputed_vcf`.
- **File? merged_qc_metrics**: Merged QC metrics TSV file, if `qc_metrics` input is defined and contains non-_null_ values.

## Glimpse2ImputationInBatches

This workflow splits the provided CRAMs into batches, performs imputation using [Glimpse2Imputation](#Glimpse2Imputation), and merges the batches back together using [Glimpse2MergeBatches](#Glimpse2MergeBatches). In addition to the inputs to imputation, this workflow takes an **Int batch_size**. This workflow can only operate on CRAM inputs, it does not support VCF inputs to imputation. 

## Appendix
### Glimpse2MergeBatches AF and INFO score recalculation
**By: Christopher Kachulis**

Glimpse outputs three values in the info field:
1) AF.  This is easily calculated for the full cohort based on the individual batch values as a weighted mean of the batch allele frequencies, $$AF_{cohort}=\frac{\sum AF_i N_i}{\sum N_i}$$

2) RAF.  This is the reference panel allele frequency.  Assuming the same reference panel was used for all batches, these are all the same, so just take the first value

3) INFO.  This is the "IMPUTE style info score".  Based on https://static-content.springer.com/esm/art%3A10.1038%2Fnrg2796/MediaObjects/41576_2010_BFnrg2796_MOESM3_ESM.pdf, this is calculated as $$1-\frac{\sum f_j - e^2_j}{2N\times AF(1-AF)}$$, or 1 if AF=0,1, with $j$ running across the N samples, $f_j = p_{j1} + 4 p_{j2}$, $e_j = p_{j1} + 2 p_{j2}$, where $p_{j1}$ is the imputed posterior of a het for sample $j$, and $p_{j2}$ is the imputed posterior of a hom var for sample $j$.  Note that the terms in the denominator are all easily calculated for a cohort based on their values for the constituent batches;  AF as described above, and N just as the sum over the batches.  The numerator we can define for batch $i$ as $C_i = \sum (f_{ij} -e_{ij}^2)$ and note that, for a whole cohort, we simply have $C_{cohort} = \sum C_i$.  We then note that, for a batch i, we have $$I_i =1 - \frac{C_i}{2N_i\times AF_i(1-AF_i)}$$, so we can solve for $C_i$ as $$C_i = (1-I_i)*2N_i\times AF_i(1-AF_i)$$.  We can then calculate the cohort INFO score as $$I_{cohort}=1-\frac{\sum C_i}{2N_{cohort} AF_{cohort}(1-AF_{cohort})}$$ which becomes, $$I_{cohort}=1-\frac{\sum (1-I_i)*2N_i\times AF_i(1-AF_i)}{2\sum N_i \times \frac{\sum AF_i N_i}{\sum N_i} (1-\frac{\sum AF_i N_i}{\sum N_i})}$$
