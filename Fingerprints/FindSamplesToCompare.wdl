version 1.0
import "../BenchmarkVCFs/BenchmarkVCFs.wdl" as Benchmark

workflow FindSamplesToCompare {
    
    input {
        Array[File] input_callset
        Array[File] ground_truth_files
        Array[File] ground_truth_indexes
        Array[File] ground_truth_intervals
        Array[String] truth_labels
        Array[File] annotation_intervals
  
        File ref_fasta
        File ref_fasta_index
        File ref_fasta_dict
        File ref_fasta_sdf

        File haplotype_database 
        File picard_cloud_jar

        String docker 

        File? interval_list_override
        File? runs_file_override
        String comparison_docker

        Boolean remove_symbolic_alleles=false
    }


    File interval_list = select_first([interval_list_override, "gs://concordance/interval_lists/chr9.hg38.interval_list"])
    File runs_file = select_first([runs_file_override, "gs://concordance/hg38/runs.conservative.bed"])
    

    Int VCF_disk_size = 50
    File monitoring_script="gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"

    call MakeStringMap as intervalsMap {input: keys=ground_truth_files, values=ground_truth_intervals}
    call MakeStringMap as lablesMap    {input: keys=ground_truth_files, values=truth_labels}
    call MakeStringMap as indexesMap   {input: keys=ground_truth_files, values=ground_truth_indexes}
 
    Map[File, File]   truthIntervals = intervalsMap.map
    Map[File, String] truthLabels    = lablesMap.map
    Map[File, File]   truthIndex     = indexesMap.map


    call SwitchFilterAnnotation {
        input:
            input_vcf = input_callset[0],
            INFO_TAG_OLD="VQSLOD",
            INFO_TAG_NEW="TREE_SCORE",
            output_vcf_basename = "callset_swapped_score",
            preemptible_tries = 0,
            disk_size = 2*size(input_callset, "GiB") + 20
      }

    call CrosscheckFingerprints {
         input:
           input_data = [SwitchFilterAnnotation.output_vcf],
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

            call ExtractSampleFromCallset {
                input:
                    callset=match.leftFile,
                    sample=match.leftSample,
                    basename=match.leftSample + ".extracted"
            }

            call Benchmark.Benchmark{
                input:
                     analysisRegion = "chr9",
                     evalVcf = ExtractSampleFromCallset.output_vcf,
                     evalLabel = match.leftSample,
                     evalVcfIndex = ExtractSampleFromCallset.output_vcf_index,
                     truthVcf = match.rightFile,
                     confidenceInterval = truthIntervals[match.rightFile],
                     truthLabel = match.rightSample,
                     truthVcfIndex = truthIndex[match.rightFile],
                     reference = ref_fasta,
                     refIndex = ref_fasta_index,
                     refDict = ref_fasta_dict,
                     hapMap = haplotype_database,
                     stratIntervals = [ 
                     "gs://concordance/hg38/runs.conservative.bed", 
                     "gs://concordance/hg38/LCR-hs38.bed",
                     "gs://concordance/hg38/mappability.0.bed",
                     "gs://concordance/hg38/exome.twist.bed"], ## TODO
                     stratLabels = [
                     "HMER 7+", 
                     "LCR-hs38",
                     "mappability=0",
                     "exome"], ## TODO
                     jexlVariantSelectors = ["vc.isSimpleIndel()  && vc.getIndelLengths().0<0", "vc.isSimpleIndel() && vc.getIndelLengths().0>0"],
                     variantSelectorLabels = ["deletion","insertion"],
                     referenceVersion = "1",
                     doIndelLengthStratification=false,
                     gatkTag="4.0.11.0",
                     requireMatchingGenotypes=true,
                     passingOnly=true,
                     vcfScoreField = "INFO.TREE_SCORE"
                 }

            Pair[File,File] vcf_and_index_original = zip([ExtractSampleFromCallset.output_vcf],[ExtractSampleFromCallset.output_vcf_index])[0]
            
            Float compareSize = 3*size([ExtractSampleFromCallset.output_vcf,match.rightFile,ref_fasta,ref_fasta_sdf], "GiB")

            if (remove_symbolic_alleles){
                call FilterSymbolicAlleles{
                    input:
                        monitoring_script=monitoring_script,
                        input_vcf=ExtractSampleFromCallset.output_vcf,
                        input_vcf_index=ExtractSampleFromCallset.output_vcf_index,
                        output_basename=match.leftSample + ".noSymbolicAlleles",
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        ref_fasta_dict = ref_fasta_dict,
                        disk_size = round(compareSize),
                        preemptible_tries = 3,
                        no_address = true,
                }

                Pair[File,File] vcf_and_index_symbolic_removed = zip([FilterSymbolicAlleles.output_vcf],[FilterSymbolicAlleles.output_vcf_index])[0]
            }


            Pair[File,File] vcf_and_index_to_compare = select_first([vcf_and_index_symbolic_removed,vcf_and_index_original])

            call CompareToGroundTruth {
                input:
                    monitoring_script=monitoring_script,
                    left_sample_name = match.leftSample,
                    right_sample_name = match.rightSample,
                    input_vcf = vcf_and_index_to_compare.left,
                    input_vcf_index = vcf_and_index_to_compare.right,
                    input_vcf_name = match.leftSample,
                    truth_vcf = match.rightFile,
                    truth_vcf_index = truthIndex[match.rightFile],
                    truth_high_confidence = truthIntervals[match.rightFile],
                    interval_list=interval_list,
                    runs_file=runs_file,
                    ref_fasta=ref_fasta,
                    ref_fasta_sdf=ref_fasta_sdf,
                    annotation_intervals=annotation_intervals,
                    preemptible_tries=3,
                    disk_size=compareSize+20,
                    docker=comparison_docker,
                    no_address=true
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
        File output_vcf_index = "~{basename}.vcf.gz.tbi"
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
    File truth_vcf 
    File truth_vcf_index 
    File truth_high_confidence

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
            --gtr_vcf ~{truth_vcf} \
            --cmp_intervals ~{interval_list} \
            --highconf_intervals ~{truth_high_confidence} \
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


task FilterSymbolicAlleles { 
  input {
    File monitoring_script
    File input_vcf
    File input_vcf_index
    String output_basename
    Int disk_size
    File ref_fasta
    File ref_fasta_index
    File ref_fasta_dict
    Int preemptible_tries
    Boolean no_address
  }
  command <<<
    bash ~{monitoring_script} > monitoring.log &
   
    set -e

    gatk --java-options "-Xmx10g" LeftAlignAndTrimVariants \
        -V  ~{input_vcf} \
        -O ~{output_basename}.tmp1.vcf.gz \
        --reference ~{ref_fasta} \
        --split-multi-allelics  

    rm ~{input_vcf}

    gatk --java-options "-Xmx10g"  SelectVariants \
          -V ~{output_basename}.tmp1.vcf.gz \
          -O ~{output_basename}.tmp2.vcf.gz \
          --remove-unused-alternates

    rm ~{output_basename}.tmp1.vcf.gz

    gatk --java-options "-Xmx10g" SelectVariants \
          -V ~{output_basename}.tmp2.vcf.gz \
          -O ~{output_basename}.vcf.gz \
          --exclude-non-variants \
          --select-type-to-exclude SYMBOLIC
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "12 GB"
    cpu: "1"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: "broadinstitute/gatk:4.1.4.1" 
    noAddress: no_address
    maxRetries: 1
  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
    File monitoring_log = "monitoring.log"
  }
}



task SwitchFilterAnnotation {
  input {
    File input_vcf
    String INFO_TAG_OLD
    String INFO_TAG_NEW
    String output_vcf_basename
    Int preemptible_tries
    Float disk_size
  }

  command <<<
    set -xe

    # assumes gziped file 
    zgrep -m1 "##INFO=<ID=~{INFO_TAG_OLD}," ~{input_vcf} | sed 's/~{INFO_TAG_OLD}/~{INFO_TAG_NEW}/' > new_header_line.txt
    

    mkdir links
    ln -s ~{input_vcf} links/
    local_vcf=links/$(basename ~{input_vcf})

    bcftools index -t ${local_vcf}
        
    bcftools annotate -a ${local_vcf} -h new_header_line.txt -c "+~{INFO_TAG_NEW}:=~{INFO_TAG_OLD}"  ${local_vcf} -O z -o ~{output_vcf_basename}.scored.vcf.gz
    bcftools index -t ~{output_vcf_basename}.scored.vcf.gz

  >>>
  output {
    File output_vcf="~{output_vcf_basename}.scored.vcf.gz"
    File output_vcf_index="~{output_vcf_basename}.scored.vcf.gz.tbi"
  }
  runtime {
    preemptible: preemptible_tries
    memory: "16 GB"
    disks: "local-disk " + ceil(disk_size + 20) + " SSD"
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    noAddress: true
  }
}
