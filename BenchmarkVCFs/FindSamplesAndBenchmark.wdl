version 1.0
import "SimpleBenchmark.wdl" as Benchmark
import "BenchmarkVCFs.wdl" as Tasks

workflow FindSamplesAndBenchmark {

    input {
        Array[File] input_callset
        Array[String] input_callset_labels
        Array[File] ground_truth_files
        Array[File] ground_truth_indexes
        Array[File] ground_truth_intervals
        Array[String] truth_labels


        #####################
        ## Optional Inputs ##
        #####################
        File? gatkJarForAnnotation
        Array[String]? annotationNames = []

        File ref_fasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
        File ref_fasta_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        File ref_fasta_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
        File haplotype_database = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt"
        String referenceVersion = "hg38"

        String vcf_score_field = "INFO.TREE_SCORE"

        String gatkTag = "4.2.3.0"
        File picard_cloud_jar = "gs://broad-dsde-methods/picard/picardcloud-2.26.6.jar"
        String picardDocker = "us.gcr.io/broad-gatk/gatk:4.2.3.0"
        String bcftoolsDocker = "us.gcr.io/broad-dsde-methods/imputation_bcftools_vcftools_docker:v1.0.0"
        String pythonDocker = "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"

        String? analysis_region

        Boolean remove_symbolic_alleles = false
        Boolean passingOnly = true
        Boolean doIndelLengthStratification = false
        Boolean requireMatchingGenotypes = true

        Array[File] stratIntervals = []
        Array[String] stratLabels = []

        Array[String] jexlVariantSelectors = ["vc.isSimpleIndel()  && vc.getIndelLengths().0<0", "vc.isSimpleIndel() && vc.getIndelLengths().0>0"]
        Array[String] variantSelectorLabels = ["deletion","insertion"]

        Array[String]? experiment_label
        Array[String]? extra_column
        String? extra_column_label

        # Input for monitoring_script can be found here: https://github.com/broadinstitute/palantir-workflows/blob/main/Scripts/monitoring/cromwell_monitoring_script.sh.
        # It must be copied to a google bucket and then the google bucket path can be used as the input for monitoring_script.
        File monitoring_script = "gs://broad-dsde-methods-hydro-gen-truth-data-public/scripts/cromwell_monitoring_script.sh"

        String? dummyInputForTerraCallCaching
    }

    parameter_meta {
        input_callset_labels: {description:"Labels for each input callset that is being evaluated. This will include this label in the output summary table under the 'Name' column."}
    }

    Int VCF_disk_size = ceil(size(input_callset, "GiB")) + 10

    String gatkDocker = "us.gcr.io/broad-gatk/gatk:" + gatkTag

    call MakeStringMap as intervalsMap {input: keys=ground_truth_files, values=ground_truth_intervals}
    call MakeStringMap as labelsMap    {input: keys=ground_truth_files, values=truth_labels}
    call MakeStringMap as indexesMap   {input: keys=ground_truth_files, values=ground_truth_indexes}
    call MakeStringMap as evalLabelsMap {input: keys=input_callset,     values=input_callset_labels}
    if (defined(experiment_label)) {call MakeStringMap as experimentLabelMap {input: keys=input_callset, values=select_first([experiment_label])}}
    if (defined(extra_column))     {call MakeStringMap as extraColumnMap     {input: keys=input_callset, values=select_first([extra_column])}}

    Map[File, File]   truthIntervals = intervalsMap.map
    Map[File, String] truthLabels    = labelsMap.map
    Map[File, File]   truthIndex     = indexesMap.map
    Map[File, String] evalLabels     = evalLabelsMap.map
    Map[File, String]? experimentLabels = experimentLabelMap.map
    Map[File, String]? extraColumns     = extraColumnMap.map

    call CrosscheckFingerprints {
        input:
            input_data = input_callset,
            metrics_basename = "crosscheck",
            ground_truth_files = ground_truth_files,
            ground_truth_indexes = ground_truth_indexes,
            haplotype_database = haplotype_database,
            disk_size = VCF_disk_size,
            preemptible_tries = 3,
            picardDocker = picardDocker,
            monitoring_script = monitoring_script,
            picard_jar = picard_cloud_jar,
    }

    call PickMatches {
        input:
            crosscheck_results = CrosscheckFingerprints.crosscheck,
            ground_truth_global = CrosscheckFingerprints.ground_truth_global,
            ground_truth_local = CrosscheckFingerprints.ground_truth_local,
            pythonDocker = pythonDocker
    }


    scatter(interval_and_label in zip(stratIntervals, stratLabels)){
        call Tasks.ConvertIntervals as ConvertIntervals {
            input:
                inputIntervals = interval_and_label.left,
                refDict = ref_fasta_dict,
                gatkTag = gatkTag,
                dummyInputForTerraCallCaching = dummyInputForTerraCallCaching
        }

        call InvertIntervalList {
            input:
                interval_list = ConvertIntervals.intervalList,
                picardDocker = picardDocker,
                picard_jar = picard_cloud_jar,
                disk_size = VCF_disk_size,
                dummyInputForTerraCallCaching = dummyInputForTerraCallCaching
        }
        String notLabel="NOT_" + interval_and_label.right
    }

    Array[File] allStratIntervals = flatten([stratIntervals,InvertIntervalList.output_interval])
    Array[String] allStratLabels = flatten([stratLabels,notLabel])

    scatter(matchArray in PickMatches.matches) {
        Match match = object{
                          leftFile: matchArray[0],
                          leftSample: matchArray[1],
                          rightFile: matchArray[2],
                          rightSample: matchArray[3]
                      }

        if defined(experimentLabels) {String? thisExperimentLabel = select_first([experimentLabels])[match.rightFile]}
        if defined(extraColumns)     {String? thisExtraColumn     = select_first([extraColumns])[match.rightFile]}

        call ExtractSampleFromCallset {
            input:
                callset = match.leftFile,
                sample = match.leftSample,
                basename = match.leftSample + ".extracted",
                bcftoolsDocker = bcftoolsDocker
        }

        call Benchmark.SimpleBenchmark as BenchmarkVCF{
            input:
                base_vcf = match.rightFile,
                base_vcf_index = truthIndex[match.rightFile],
                base_output_sample_name = truthLabels[match.rightFile],
                base_vcf_sample_name = match.rightSample,
                call_vcf = ExtractSampleFromCallset.output_vcf,
                call_vcf_index = ExtractSampleFromCallset.output_vcf_index,
                call_output_sample_name = match.leftSample,
                call_vcf_sample_name = match.leftSample,
                ref_fasta = ref_fasta,
                ref_index = ref_fasta_index,
                ref_dict = ref_fasta_dict,
                strat_intervals = allStratIntervals,
                strat_labels = allStratLabels,
                bcf_selectors = jexlVariantSelectors,
                bcf_labels = variantSelectorLabels,
                evaluation_bed = truthIntervals[match.rightFile],
                experiment_label = thisExperimentLabel,
                extra_column = thisExtraColumn,
                extra_column_label = extra_column_label,
                preemptible = 1
        }

        Pair[File,File] vcf_and_index_original = zip([ExtractSampleFromCallset.output_vcf],[ExtractSampleFromCallset.output_vcf_index])[0]

        Float compareSize = 3*size([ExtractSampleFromCallset.output_vcf,match.rightFile,ref_fasta], "GiB")

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
                    gatkDocker = gatkDocker
            }

            Pair[File,File] vcf_and_index_symbolic_removed = zip([FilterSymbolicAlleles.output_vcf],[FilterSymbolicAlleles.output_vcf_index])[0]
        }


        Pair[File,File] vcf_and_index_to_compare = select_first([vcf_and_index_symbolic_removed,vcf_and_index_original])

    }
    call Tasks.CombineSummaries as CombineSummariesROC{
        input:
            summaries = BenchmarkVCF.combined_ROC,
            delimeter = "\t",
            quotes_in_output = false,
            output_filename = "Full_ROC.tsv"
    }
    call Tasks.CombineSummaries as CombineSummariesSimple{
        input:
            summaries = BenchmarkVCF.simple_summary,
            delimeter = "\t",
            quotes_in_output = false,
            output_filename = "SimpleSummary.tsv"
    }
    call Tasks.CombineSummaries as CombineSummariesIDD{
        input:
            summaries = BenchmarkVCF.combined_IDD,
            delimeter = "\t",
            quotes_in_output = false,
            output_filename = "Full_IDD.tsv"
    }
    call Tasks.CombineSummaries as CombineSummariesST{
        input:
            summaries = BenchmarkVCF.combined_ST,
            delimeter = "\t",
            quotes_in_output = false,
            output_filename = "Full_ST.tsv"
    }

    output {
        File benchmark_summary = CombineSummariesSimple.summaryOut
        File benchmark_ROC = CombineSummariesROC.summaryOut
        File benchmark_IDD = CombineSummariesIDD.summaryOut
        File benchmark_ST = CombineSummariesST.summaryOut
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
        Array[File] ground_truth_indexes
        File haplotype_database
        Int disk_size
        Int preemptible_tries
        String picardDocker
        File picard_jar
    }
    parameter_meta {
        input_data: {
                        localization_optional: true
                    }

        ground_truth_files: {
                                localization_optional: false
                            }
    }

    String tsv_out="~{metrics_basename}.fingerprints_metrics"
    Array[String] ground_truth_strings = ground_truth_files

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

        # Put local Cromwell paths into one file, and global bucket paths into another
        echo "LOCAL" > ground_truth_files.csv
        echo "~{sep="\n" ground_truth_files}" >> ground_truth_files.csv
        echo "GLOBAL" > ground_truth_strings.csv
        echo "~{sep="\n" ground_truth_strings}" >> ground_truth_strings.csv
    >>>
    output {
        File crosscheck = tsv_out
        File ground_truth_local = "ground_truth_files.csv"
        File ground_truth_global = "ground_truth_strings.csv"
    }

    runtime {
        preemptible: preemptible_tries
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: picardDocker
        noAddress: false
        maxRetries: 1
        continueOnReturnCode: [0,5]
    }
}

task PickMatches {
    input {
        File crosscheck_results
        File ground_truth_local
        File ground_truth_global
        String pythonDocker
    }

    command <<<
        python <<CODE
        import pandas as pd

        # Import crosscheck results into DataFrame
        crosscheck = pd.read_csv("~{crosscheck_results}", sep='\t', comment='#', skip_blank_lines=True)

        # Find the pairs which were (expected/unexpected) matches
        match_list = crosscheck[-crosscheck['RESULT'].str.contains('MISMATCH')]

        # Fix file path back to cloud rather than localized version
        truth_files_local = pd.read_csv("~{ground_truth_local}")
        truth_files = pd.read_csv("~{ground_truth_global}")
        truth_files['LOCAL'] = "file://" + truth_files_local   # Fixes some formatting issue from the crosscheck output
        local_to_global = dict(zip(truth_files.LOCAL, truth_files.GLOBAL))
        match_list['RIGHT_FILE'] = match_list['RIGHT_FILE'].map(local_to_global)

        # Write out the left/right files/samples from the list of matches
        match_list[['LEFT_FILE', 'LEFT_SAMPLE', 'RIGHT_FILE', 'RIGHT_SAMPLE']].to_csv("matches.tsv", sep='\t', index=False, header=False)

        CODE
    >>>

    output {
        Array[Array[String]] matches = read_tsv("matches.tsv")
    }

    runtime {
        docker: pythonDocker
    }
}



task ExtractSampleFromCallset {
    input {
        File callset
        String sample
        String basename
        String bcftoolsDocker
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
        String gatkDocker
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
        docker: gatkDocker
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
        String picardDocker
        File picard_jar
        Int disk_size
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
        disks: "local-disk " + disk_size + " HDD"
        docker: picardDocker
    }
}
