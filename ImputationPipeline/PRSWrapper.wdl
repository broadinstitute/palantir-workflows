version 1.0
import "PRSTasks.wdl" as PRSTasks
import "ScoringPart.wdl" as Score

workflow PRSWrapper {
  input {
    Array[String] condition_names
    Array[Boolean] use_condition
    Array[Boolean] use_ancestry_correction
    Array[File] weights_files

    File vcf
    File population_vcf
  }

  if (length(condition_names) != length(use_condition) || length(condition_names) != length(use_ancestry_correction) || length(condition_names) != length(weights_files)) {
    call PRSTasks.ErrorWithMessage {
      input:
        message = "conditions_names, use_condition, use_ancestry_correction, and weights_files must all be arrays of the same length"
    }
  }

  scatter(i in range(length(condition_names))) {

  }
}
