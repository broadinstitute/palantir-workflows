version 1.0

# Legacy tasks from old BenchmarkSVs pipeline, preserved for anyone looking to experiment with Wittyer

# Base for this task was originally written by Yueyao Gao here:
# https://github.com/broadinstitute/TAG-public/blob/yg_subset_callset/BenchmarkCNV/BenchmarkCNV/BenchmarkCNV.wdl
task WittyerEval {
    input {
        File truth_vcf
        Array[String]? truth_sample_names

        File eval_vcf
        Array[String]? query_sample_names

        File? bedfile

        File? wittyer_config
        String wittyer_evaluation_mode = "CrossTypeAndSimpleCounting"

        RuntimeAttributes runtimeAttributes = {"disk_size": 100, "cpu": 4, "memory": 8}
    }

    command <<<
        set -x

        ### Create sample pair list
        TRUTH_SAMPLES=(~{sep=' ' truth_sample_names})
        QUERY_SAMPLES=(~{sep=' ' query_sample_names})

        # Zip names together
        unset ZIPPED_LIST
        for (( i=0; i<${#TRUTH_SAMPLES[*]}; ++i)); do
            ZIPPED_LIST+=( "${TRUTH_SAMPLES[$i]}:${QUERY_SAMPLES[$i]}" )
        done

        PAIRS_STRING=$(IFS=, ; echo "${ZIPPED_LIST[*]}")

        # Create flag
        if (( ${#ZIPPED_LIST[@]} )); then
            SP_FLAG="-sp ${PAIRS_STRING}"
        else
            SP_FLAG=""
        fi

        # Run Benchmarking tool wittyer
        /opt/Wittyer/Wittyer \
            -i ~{eval_vcf} \
            -t ~{truth_vcf} \
            -em ~{wittyer_evaluation_mode} \
            ~{"--configFile " + wittyer_config} \
            ~{'--includeBed ' + bedfile} \
            -o wittyer_output \
            $SP_FLAG

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/wittyer:v1.0"
        bootDiskSizeGb: 12
        memory: runtimeAttributes.memory + " GB"
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        cpu: runtimeAttributes.cpu
        preemptible: 1
    }

    output {
        File wittyer_stats = "wittyer_output/Wittyer.Stats.json"
        File wittyer_config = "wittyer_output/Wittyer.ConfigFileUsed.json"
        Array[File] wittyer_annotated_vcf = glob("wittyer_output/Wittyer.*.Vs.*.vcf.gz")
        Array[File] wittyer_annotated_vcf_index = glob("wittyer_output/Wittyer.*.Vs.*.vcf.gz.tbi")
    }
}

task AddIntervalOverlapStatsWittyer {
    input {
        File input_vcf

        String? base_sample_name
        String? comp_sample_name
        String? experiment

        Array[File] bed_regions
        Array[String] bed_labels
        Int breakpoint_padding = 20

        File ref_fai

        # Runtime parameters
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(input_vcf, "GB")) + 100, "cpu": 4, "memory": 16}
    }

    command <<<
        set -xueo pipefail

        # Extract useful columns from Wittyer VCF
        # Use POS0/END0 for conherence with bed coordinates later
        bcftools query -i 'INFO/SVTYPE!="."' -f'%CHROM\t%POS0\t%END0\t%SVLEN\t%SVTYPE\t%FILTER\t%WHERE\t%WIN[\t%GT\t%WIT\t%WHY]\n' ~{input_vcf} | bedtools sort -i - > wittyer_vcf_labels.tsv

        # Format table in Python
        python3 << CODE
        import pandas as pd

        df = pd.read_csv('wittyer_vcf_labels.tsv', sep='\t', header=None)
        columns = ['CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'FILTER', 'WHERE', 'WIN', 'GT', 'WIT', 'WHY']

        no_gt = df[(df[8] == '.') & (df[11] == '.')].drop(columns=[8, 11])
        truth_df = df[df[8] != '.'].drop(columns=[11, 12, 13])
        query_df = df[df[11] != '.'].drop(columns=[8, 9, 10])

        no_gt.columns = ['CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'FILTER', 'WHERE', 'WIN', 'TWIT', 'TWHY', 'QWIT', 'QWHY']
        truth_df.columns = columns
        query_df.columns = columns

        truth_df['TruthSample'] = "~{base_sample_name}"
        truth_df['QuerySample'] = "~{comp_sample_name}"
        query_df['TruthSample'] = "~{base_sample_name}"
        query_df['QuerySample'] = "~{comp_sample_name}"

        truth_df['Experiment'] = "~{experiment}"
        query_df['Experiment'] = "~{experiment}"
        no_gt['Experiment'] = "~{experiment}"

        no_gt.to_csv('no_gt.tsv', sep='\t', index=False, header=False)
        truth_df.to_csv('truth_table.tsv', sep='\t', index=False, header=False)
        query_df.to_csv('query_table.tsv', sep='\t', index=False, header=False)

        CODE

        # Clean ref_fai into genome file for bedtools
        cut -f1,2 ~{ref_fai} > ref.genome

        # Collect interval bed overlap stats
        generate_interval_stats () {
            CURRENT_LABEL=$1
            INTERVAL_FILE=$2
            INPUT_LABEL=$3
            VCF_FILE=$4

            echo -e "${CURRENT_LABEL}-count\t${CURRENT_LABEL}-overlap" > header.txt
            # WARNING: annotate does not preserve order of input bed regions!
            # Also: bedtools annotate complains with this input if not cut to first 3 fields
            cut -f1-3 $VCF_FILE > vcf.bed
            bedtools annotate -both -i vcf.bed -files $INTERVAL_FILE | bedtools sort -i - | rev | cut -f1-2 | rev > "${CURRENT_LABEL}-${INPUT_LABEL}-preheader.bed"
            cat header.txt "${CURRENT_LABEL}-${INPUT_LABEL}-preheader.bed" > "${CURRENT_LABEL}-${INPUT_LABEL}-annotated.bed"

            # Also add in breakpoint overlap stats using POS and END values
            echo -e "${CURRENT_LABEL}-LBEND-count\t${CURRENT_LABEL}-LBEND-overlap" > pos-header.txt
            echo -e "${CURRENT_LABEL}-RBEND-count\t${CURRENT_LABEL}-RBEND-overlap" > end-header.txt
            echo "Printing vcf.bed head..."
            head -n 10 vcf.bed
            awk '{OFS="\t"; print $1, $2, $2}' vcf.bed | bedtools slop -b ~{breakpoint_padding} -i - -g ref.genome > pos.bed
            awk '{OFS="\t"; print $1, $3, $3}' vcf.bed | bedtools slop -b ~{breakpoint_padding} -i - -g ref.genome > end.bed
            bedtools annotate -both -i pos.bed -files $INTERVAL_FILE | bedtools sort -i - | rev | cut -f1-2 | rev > "${CURRENT_LABEL}-${INPUT_LABEL}-pos-preheader.bed"
            bedtools annotate -both -i end.bed -files $INTERVAL_FILE | bedtools sort -i - | rev | cut -f1-2 | rev > "${CURRENT_LABEL}-${INPUT_LABEL}-end-preheader.bed"
            cat pos-header.txt "${CURRENT_LABEL}-${INPUT_LABEL}-pos-preheader.bed" > "${CURRENT_LABEL}-${INPUT_LABEL}-pos-annotated.bed"
            cat end-header.txt "${CURRENT_LABEL}-${INPUT_LABEL}-end-preheader.bed" > "${CURRENT_LABEL}-${INPUT_LABEL}-end-annotated.bed"

            # Paste together three overlap stats files
            paste "${CURRENT_LABEL}-${INPUT_LABEL}-annotated.bed" "${CURRENT_LABEL}-${INPUT_LABEL}-pos-annotated.bed" "${CURRENT_LABEL}-${INPUT_LABEL}-end-annotated.bed" > "${CURRENT_LABEL}-${INPUT_LABEL}-full-annotated.bed"
        }

        INTERVAL_FILES=(~{sep=' ' bed_regions})
        INTERVAL_LABELS=(~{sep=' ' bed_labels})

        for i in {0..~{length(bed_regions)-1}}
        do
            CURRENT_LABEL="${INTERVAL_LABELS[$i]}"
            CURRENT_FILE="${INTERVAL_FILES[$i]}"
            generate_interval_stats $CURRENT_LABEL $CURRENT_FILE "no-gt" no_gt.tsv
            generate_interval_stats $CURRENT_LABEL $CURRENT_FILE "truth-table" truth_table.tsv
            generate_interval_stats $CURRENT_LABEL $CURRENT_FILE "query-table" query_table.tsv
        done

        # Combine across interval beds
        echo -e "CHROM\tPOS\tEND\tSVLEN\tSVTYPE\tFILTER\tWHERE\tWIN\tTWIT\tTWHY\tQWIT\tQWHY\tExperiment" > no_gt-header.txt
        cat no_gt-header.txt no_gt.tsv > no_gt-with_header.tsv
        paste no_gt-with_header.tsv *-no-gt-full-annotated.bed > no-gt_intervals-final.bed

        echo -e "CHROM\tPOS\tEND\tSVLEN\tSVTYPE\tFILTER\tWHERE\tWIN\tGT\tWIT\tWHY\tTruthSample\tQuerySample\tExperiment" > header.txt
        cat header.txt truth_table.tsv > truth_table-with_header.tsv
        paste truth_table-with_header.tsv *-truth-table-full-annotated.bed > truth-table_intervals-final.bed

        cat header.txt query_table.tsv > query_table-with_header.tsv
        paste query_table-with_header.tsv *-query-table-full-annotated.bed > query-table_intervals-final.bed

        # Correct POS0/END back to original VCF POS/END coordinates
        python3 << CODE
        import pandas as pd

        for file in ['no-gt_intervals-final.bed', 'truth-table_intervals-final.bed', 'query-table_intervals-final.bed']:
            df = pd.read_csv(file, sep='\t')
            df['POS'] = df['POS'] + 1
            df['END'] = df['END'] + 1
            df['VCF'] = file.split('-')[0]    # Add truth/query labels for respective tables
            df.to_csv(file, sep='\t', index=False)

        CODE

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.0"
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        memory: runtimeAttributes.memory + " GB"
        cpu: runtimeAttributes.cpu
    }

    output {
        File no_gt_intervals = "no-gt_intervals-final.bed"
        File truth_table_intervals = "truth-table_intervals-final.bed"
        File query_table_intervals = "query-table_intervals-final.bed"
    }
}


task CleanBasicWittyerStats {
    input {
        File wittyer_stats
        String experiment = ""
    }

    command <<<
        set -xueo pipefail

        python3 << CODE
        import pandas as pd
        import json

        with open('~{wittyer_stats}') as file:
            witt = json.load(file)

        df = pd.DataFrame()

        for sample_pair in witt['PerSampleStats']:
            query_name = sample_pair['QuerySampleName']
            truth_name = sample_pair['TruthSampleName']

            for stats in sample_pair['OverallStats']:
                new_row = dict(**{'QueryName': query_name, 'TruthName': truth_name, 'SVTYPE': 'All', 'Bin': 'All'},
                              **stats)
                df = pd.concat([df, pd.DataFrame(new_row, index=[0])])

            for variants in sample_pair['DetailedStats']:
                variant_type = variants['VariantType']

                for stats in variants['OverallStats']:
                    new_row = dict(**{'QueryName': query_name, 'TruthName': truth_name, 'SVTYPE': variant_type, 'Bin': 'All'},
                              **stats)
                    df = pd.concat([df, pd.DataFrame(new_row, index=[0])])

                    for bins in variants['PerBinStats']:
                        Bin = bins['Bin']
                        for bin_stats in bins['Stats']:
                            new_row = dict(**{'QueryName': query_name, 'TruthName': truth_name, 'SVTYPE': variant_type, 'Bin': Bin},
                                        **bin_stats)
                            df = pd.concat([df, pd.DataFrame(new_row, index=[0])])

        df['Experiment'] = "~{experiment}"
        df.to_csv('wittyer_stats-cleaned.tsv', sep='\t', index=False)

        CODE

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.0"
        disks: "local-disk " + 50 + " HDD"
        memory: 4 + " GB"
        cpu: 2
    }

    output {
        File wittyer_stats_cleaned = "wittyer_stats-cleaned.tsv"
    }
}
