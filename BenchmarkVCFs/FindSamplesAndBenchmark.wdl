version 1.0
import "BenchmarkVCFs.wdl" as Benchmark

workflow FindSamplesAndBenchmark {
    
    input {
        Array[File] input_callset
        Array[String] input_callset_labels
        Array[File] ground_truth_files
        Array[File] ground_truth_indexes
        Array[File] ground_truth_intervals
        Array[String] truth_labels
        Array[File] annotation_intervals

        File? gatkJarForAnnotation
        Array[String]? annotationNames = []
  
        File ref_fasta
        File ref_fasta_index
        File ref_fasta_dict
        File ref_fasta_sdf

        File haplotype_database 
        File picard_cloud_jar

        String vcf_score_field

        String docker
        String bcftoolsDocker

        String? analysis_region

        Boolean remove_symbolic_alleles = false
        Boolean passingOnly = true
        Boolean doIndelLengthStratification = false
        Boolean requireMatchingGenotypes = true

        String referenceVersion = "1"
        String gatkTag = "4.2.3.0"

        Array[File] stratIntervals = []
        Array[String] stratLabels = []

        Array[String] jexlVariantSelectors = ["vc.isSimpleIndel()  && vc.getIndelLengths().0<0", "vc.isSimpleIndel() && vc.getIndelLengths().0>0"]
        Array[String] variantSelectorLabels = ["deletion","insertion"]
        # Input for monitoring_script can be found here: https://github.com/broadinstitute/palantir-workflows/blob/main/Scripts/monitoring/cromwell_monitoring_script.sh.
        # It must be copied to a google bucket and then the google bucket path can be used as the input for monitoring_script.
        File monitoring_script

        String? dummyInputForTerraCallCaching
    }

    parameter_meta {
      input_callset_labels: {description:"Labels for each input callset that is being evaluated. This will include this label in the output summary table under the 'Name' column."}
    }

    Int VCF_disk_size = ceil(size(input_callset, "GiB")) + 10

    call MakeStringMap as intervalsMap {input: keys=ground_truth_files, values=ground_truth_intervals}
    call MakeStringMap as lablesMap    {input: keys=ground_truth_files, values=truth_labels}
    call MakeStringMap as indexesMap   {input: keys=ground_truth_files, values=ground_truth_indexes}
    call MakeStringMap as evalLabelsMap {input: keys=input_callset,     values=input_callset_labels}

    Map[File, File]   truthIntervals = intervalsMap.map
    Map[File, String] truthLabels    = lablesMap.map
    Map[File, File]   truthIndex     = indexesMap.map
    Map[File, String] evalLabels     = evalLabelsMap.map

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


    scatter(interval_and_label in zip(stratIntervals, stratLabels)){
        call Benchmark.ConvertIntervals as ConvertIntervals {
            input:
                inputIntervals = interval_and_label.left,
                refDict = ref_fasta_dict,
                gatkTag = gatkTag,
                dummyInputForTerraCallCaching = dummyInputForTerraCallCaching
        }

        call InvertIntervalList {
            input:
                interval_list = ConvertIntervals.intervalList,
                docker = docker,
                picard_jar = picard_cloud_jar,
                dummyInputForTerraCallCaching = dummyInputForTerraCallCaching
        }
        String notLabel="NOT_" + interval_and_label.right
    }

    Array[File] allStratIntervals=flatten([stratIntervals,InvertIntervalList.output_interval])
    Array[File] allStratLabels=flatten([stratLabels,notLabel])

    scatter(matchArray in PickMatches.matches) {
            Match match = object{
                leftFile: matchArray[0],
                leftSample: matchArray[1],
                rightFile: matchArray[2],
                rightSample: matchArray[3]
            }

            call ExtractSampleFromCallset {
                input:
                    callset = match.leftFile,
                    sample = match.leftSample,
                    basename = match.leftSample + ".extracted",
                    bcftoolsDocker = bcftoolsDocker
            }

            call Benchmark.Benchmark as BenchmarkVCF{
                input:
                     analysisRegion = analysis_region,
                     evalVcf = ExtractSampleFromCallset.output_vcf,
                     evalLabel = evalLabels[match.leftFile],
                     evalVcfIndex = ExtractSampleFromCallset.output_vcf_index,
                     truthVcf = match.rightFile,
                     confidenceInterval = truthIntervals[match.rightFile],
                     truthLabel = match.rightSample,
                     truthVcfIndex = truthIndex[match.rightFile],
                     reference = ref_fasta,
                     refIndex = ref_fasta_index,
                     refDict = ref_fasta_dict,
                     hapMap = haplotype_database,
                     stratIntervals = allStratIntervals,
                     stratLabels = allStratLabels, 
                     jexlVariantSelectors = jexlVariantSelectors,
                     variantSelectorLabels = variantSelectorLabels,
                     referenceVersion = referenceVersion,
                     doIndelLengthStratification = doIndelLengthStratification,
                     gatkTag = gatkTag,
                     requireMatchingGenotypes = requireMatchingGenotypes,
                     passingOnly = passingOnly,
                     vcfScoreField = vcf_score_field,
                     gatkJarForAnnotation = gatkJarForAnnotation,
                     annotationNames = annotationNames
                 }

            Pair[File,File] vcf_and_index_original = zip([ExtractSampleFromCallset.output_vcf],[ExtractSampleFromCallset.output_vcf_index])[0]
            
            Float compareSize = 3*size([ExtractSampleFromCallset.output_vcf,match.rightFile,ref_fasta,ref_fasta_sdf], "GiB")

            if (remove_symbolic_alleles){
                call FilterSymbolicAlleles{
                    input:
                        monitoring_script = monitoring_script,
                        input_vcf = ExtractSampleFromCallset.output_vcf,
                        input_vcf_index = ExtractSampleFromCallset.output_vcf_index,
                        output_basename = match.leftSample + ".noSymbolicAlleles",
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        ref_fasta_dict = ref_fasta_dict,
                        disk_size = round(compareSize),
                        preemptible_tries = 3,
                        no_address = true,
                        gatkTag = gatkTag
                }

                Pair[File,File] vcf_and_index_symbolic_removed = zip([FilterSymbolicAlleles.output_vcf],[FilterSymbolicAlleles.output_vcf_index])[0]
            }


            Pair[File,File] vcf_and_index_to_compare = select_first([vcf_and_index_symbolic_removed,vcf_and_index_original])

    }
    call Benchmark.CombineSummaries as CombineSummaries{
        input:
            summaries = select_all(BenchmarkVCF.summary),
            preemptible = 1
    }

   output {
    File benchmark_vcf_summary = CombineSummaries.summaryOut
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

    
    java -Xmx7g -jar ~{picard_jar} CrosscheckFingerprints \
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
        String bcftoolsDocker = "us.gcr.io/broad-dsde-methods/imputation_bcftools_vcftools_docker:v1.0.0"
    }
    Int disk_size = ceil(size(callset, "GB") + 30)
    File sampleNamesToExtract = write_lines([sample])
    command <<<
        set -xe
        mkdir out_dir
        bcftools +split ~{callset} -Oz -o out_dir -S ~{sampleNamesToExtract} -i'GT="alt"'
        mv out_dir/~{sample}.vcf.gz ~{basename}.vcf.gz
        bcftools index -t ~{basename}.vcf.gz
    >>>
    output {
        File output_vcf = "~{basename}.vcf.gz"
        File output_vcf_index = "~{basename}.vcf.gz.tbi"
    }
    runtime{
        disks: "local-disk " + disk_size + " LOCAL"
        cpu: 1
        memory: 5 + " GB"
        docker: bcftoolsDocker
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
    String gatkTag
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
    docker: "us.gcr.io/broad-gatk/gatk:"+gatkTag
    noAddress: no_address
    maxRetries: 1
  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
    File monitoring_log = "monitoring.log"
  }
}

task InvertIntervalList{
    input {
        File interval_list
        String docker
        File picard_jar
        String? dummyInputForTerraCallCaching
    }

    command <<<
        java -Xmx5g -jar ~{picard_jar} IntervalListTools \
            I=~{interval_list} \
            O=~{"NOT" + basename(interval_list)} \
            INVERT=true
    >>>
    output {
        File output_interval = "NOT" + basename(interval_list)
    }
    runtime {
        memory: "6GB"
        disks: "local-disk 20 SSD"
        docker: docker
    }
}
