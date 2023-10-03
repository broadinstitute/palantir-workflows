version 1.0
import "PRSWrapper.wdl" as Wrapper
import "Structs.wdl"

workflow SampleSetPRSWrapper {
  input {
    Array[PRSWrapperConditionResource] condition_resources
    File ckd_risk_alleles
    Float z_score_reportable_range

    File vcf
    Array[String] sample_ids
    String lab_batch_id
    Boolean is_control_sample_in
    Boolean redoPCA = false

    File population_loadings
    File population_meansd
    File population_pcs
    File pruning_sites_for_pca # and the sites used for PCA
    Int mem_extract
    Int mem_vcf_to_plink
  }
  scatter(sample_id in sample_ids) {
    call Wrapper.PRSWrapper {
      input:
        condition_resources = condition_resources,
        ckd_risk_alleles = ckd_risk_alleles,
        z_score_reportable_range = z_score_reportable_range,
        vcf = vcf,
        sample_id = sample_id,
        lab_batch_id = lab_batch_id,
        is_control_sample_in = is_control_sample_in,
        redoPCA = redoPCA,
        population_loadings = population_loadings,
        population_meansd = population_meansd,
        population_pcs = population_pcs,
        pruning_sites_for_pca = pruning_sites_for_pca,
        mem_extract = mem_extract,
        mem_vcf_to_plink = mem_vcf_to_plink
    }
  }

  output {
    Array[File] results = PRSWrapper.results
    Array[File] pcs = PRSWrapper.pcs
    Array[File] missing_sites_shifts = PRSWrapper.missing_sites_shifts
  }
}
