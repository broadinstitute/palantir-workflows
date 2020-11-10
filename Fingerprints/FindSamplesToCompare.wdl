version 1.0

workflow FindSamplesToCompare {
    
    input {
        Array[File] input_callset
        Array[File] ground_truth_files
        Array[File] ground_truth_intervals
        Array[String] truth_labels

        File haplotype_database 
        File picard_cloud_jar

        String docker 

    }

    Int VCF_disk_size = 50
    File monitoring_script="gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"

    Array[Array[String]] truth_array = transpose([ground_truth_files, ground_truth_intervals,truth_labels])

    call MakeStringMap as intervalsMap {input: keys=ground_truth_files, values=ground_truth_intervals}
    call MakeStringMap as lablesMap {input: keys=ground_truth_files, values=truth_labels}
 
    Map[File, File] truthIntervals = intervalsMap.map
    Map[File, String] truthLabels = lablesMap.map

    call CrosscheckFingerprints {
         input:
           input_data = input_callset,
           metrics_basename = "crosscheck",
           ground_truth_files = ground_truth_files,
           haplotype_database = haplotype_database,
           disk_size = VCF_disk_size,
           preemptible_tries = 3,
           docker = docker,
           monitoring_script = monitoring_script,
           picard_jar = picard_cloud_jar,
      }

      call PickMatches {
        input:
            crosscheck_results = CrosscheckFingerprints.crosscheck
      }

      scatter(matchArray in PickMatches.matches) {
            Match match = object{
                leftFile: matchArray[0],
                leftSample: matchArray[1],
                rightFile: matchArray[2],
                rightSample: matchArray[3]
            }

            call ExtractSampleFromCallset{
                input:
                    callset=match.leftFile,
                    sample=match.leftFile,
                    basename=match.leftFile + ".extracted"
            }
      }

   output {
    File crosscheck = CrosscheckFingerprints.crosscheck
    Array[Array[String]] matches = PickMatches.matches

   }
}



struct Truth {
    File truthVcf
    File truthVcfIndex
    File confidenceIntervals
    String truthLabel
}

struct Eval {
    File evalVcf
    File evalVcfIndex
    String evalLabel
}

struct Match{
    File leftFile
    String leftSample
    File rightFile
    String rightSample
}

task CrosscheckFingerprints {
  input {
    File monitoring_script
    Array[File] input_data
    String metrics_basename
    Array[File] ground_truth_files
    File haplotype_database
    Int disk_size
    Int preemptible_tries
    String docker
    File picard_jar
  }
  parameter_meta {
    input_data: {
      localization_optional: true
    }
   
    ground_truth_files: {
      localization_optional: true
    }
  }  
  
  String tsv_out="~{metrics_basename}.fingerprints_metrics"

  command <<<
    bash ~{monitoring_script} > /cromwell_root/monitoring.log &

    
    java -jar ~{picard_jar} CrosscheckFingerprints \
      INPUT=~{sep=" INPUT=" input_data} \
      SECOND_INPUT=~{sep=" SECOND_INPUT=" ground_truth_files} \
      HAPLOTYPE_MAP=~{haplotype_database} \
      LOD_THRESHOLD=-5 \
      OUTPUT=~{tsv_out} \
      CROSSCHECK_BY=FILE \
      CROSSCHECK_MODE=CHECK_ALL_OTHERS \
      TEST_INPUT_READABILITY=false \
      EXIT_CODE_WHEN_MISMATCH=5
      
   >>>
   output {
    File crosscheck = tsv_out
   }

   runtime {
    preemptible: preemptible_tries
    memory: "8 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    noAddress: false
    maxRetries: 1
    continueOnReturnCode: [0,5]
  }
}

task PickMatches {
    input {
        File crosscheck_results
    }

    command <<<

        cat <<- 'AWK' > prog.awk
        BEGIN{
            OFS="	";
            parsedhead=0;
        }
        /^#/{next} 
        /^$/{next} 
        parsedhead==0 {
            parsedhead=1;
            for (i=1; i<=NF; i++) {
                f[$i] = i
            };
            next;
            } 
        $(f["RESULT"])!~/MISMATCH/ && $(f["RESULT"])~/MATCH/ {
            print $(f["LEFT_FILE"]), $(f["LEFT_SAMPLE"]), $(f["RIGHT_FILE"]), $(f["RIGHT_SAMPLE"])
            } 
        AWK


        awk -F'\t' -f prog.awk ~{crosscheck_results} > matches.tsv
    >>>
    output {
        Array[Array[String]] matches = read_tsv("matches.tsv")
    }
    runtime {
        docker: "ubuntu"
    }
}

task ExtractSampleFromCallset { 
    input {
        File callset
        String sample
        String basename
    }
    command <<<

        set -xe

        # Silly GATK needs an indexed file....
        gatk --java-options "-Xmx4g"  \
            IndexFeatureFile \
            -I ~{callset} 

        gatk --java-options "-Xmx4g"  \
            SelectVariants \
            -V ~{callset} \
            -sn ~{sample} \
            -O ~{basename}.vcf.gz
    >>>
    output {
        File output_vcf = "~{basename}.vcf.gz"
        File? output_vcf_index = "~{basename}.vcf.gz.tbx"
    } 
    runtime{
        disks: "local-disk " + 40 + " LOCAL"
        cpu: 1
        memory: 5 + " GB"
        docker: "broadinstitute/gatk:4.1.4.1"
    }   
}


task CompareToGroundTruth {
  input {
    File monitoring_script
#    File fingerprints_csv
    String left_sample_name
    String right_sample_name
    File input_vcf
    File input_vcf_index
    String input_vcf_name
    File interval_list
    File runs_file
    File ref_fasta
    File ref_fasta_sdf
    Array[File] annotation_intervals
    Int preemptible_tries
    Float disk_size
    String docker
    Boolean no_address
  }
  command <<<
    bash ~{monitoring_script} > /cromwell_root/monitoring.log &

    
    source ~/.bashrc
    conda activate genomics.py3

    python -m tarfile -e ~{ref_fasta_sdf} ~{ref_fasta}.sdf

    cd /VariantCalling/src/python/pipelines
    python run_comparison_pipeline.py \
            --n_parts 0 \
            --hpol_filter_length_dist 12 10 \
            --input_prefix $(echo "~{input_vcf}" | sed 's/\(.vcf.gz\|.vcf\)$//') \
            --output_file /cromwell_root/~{input_vcf_name}.comp.h5 \
            --gtr_vcf $(find $ground_truths_dirname | grep -E .vcf.gz$) \
            --cmp_intervals ~{interval_list} \
            --highconf_intervals $(find $ground_truths_dirname | grep -E .bed$) \
            --runs_intervals ~{runs_file} \
            --reference ~{ref_fasta} \
            --call_sample_name ~{left_sample_name} \
            --ignore_filter_status \
            --truth_sample_name ~{right_sample_name}\
            --annotate_intervals ~{sep=" --annotate_intervals " annotation_intervals} \
            --output_suffix ''
    >>>
  runtime {
    preemptible: preemptible_tries
    memory: "16 GB"
    disks: "local-disk " + ceil(disk_size) + " SSD"
    docker: docker
    noAddress: no_address
  }
  output {
    File compare_h5 = "/cromwell_root/~{input_vcf_name}.comp.h5"
    File monitoring_log = "/cromwell_root/monitoring.log"
    File snp_fp_bed = "/cromwell_root/~{input_vcf_name}.comp_snp_fp.bed"
    File snp_fn_bed = "/cromwell_root/~{input_vcf_name}.comp_snp_fn.bed"
    File hmer_fp_1_3_bed = "/cromwell_root/~{input_vcf_name}.comp_hmer_fp_1_3.bed"
    File hmer_fn_1_3_bed = "/cromwell_root/~{input_vcf_name}.comp_hmer_fn_1_3.bed"
    File hmer_fp_4_7_bed = "/cromwell_root/~{input_vcf_name}.comp_hmer_fp_4_7.bed"
    File hmer_fn_4_7_bed = "/cromwell_root/~{input_vcf_name}.comp_hmer_fn_4_7.bed"
    File hmer_fp_8_end_bed = "/cromwell_root/~{input_vcf_name}.comp_hmer_fp_8_end.bed"
    File hmer_fn_8_end_bed = "/cromwell_root/~{input_vcf_name}.comp_hmer_fn_8_end.bed"
    File non_hmer_fp_bed = "/cromwell_root/~{input_vcf_name}.comp_non_hmer_fp.bed"
    File non_hmer_fn_bed = "/cromwell_root/~{input_vcf_name}.comp_non_hmer_fn.bed"
    File genotyping_errors_fp_bed = "/cromwell_root/~{input_vcf_name}.comp_genotyping_errors_fp.bed"
    File genotyping_errors_fn_bed = "/cromwell_root/~{input_vcf_name}.comp_genotyping_errors_fn.bed"
  }
}

# only works with simple maps (two columns)
task MakeStringMap {
    input {
        Array[String] keys
        Array[String] values
    }

    String results_path = "results.tsv"
    command <<<
       cat ~{write_tsv(transpose([keys,values]))} > ~{results_path} 
    >>>
    runtime {
        docker: "python:3"
    }
    output {
        Map[String, String] map = read_map(results_path)
    }
}




