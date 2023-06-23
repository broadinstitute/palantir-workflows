version 1.0

workflow BenchmarkSVs {
    input {
        File? base_vcf
        File? base_vcf_index
        Array[String] base_sample_names

        File comp_vcf
        File comp_vcf_index
        Array[String] comp_sample_names

        String? experiment

        File ref_dict

        Boolean split_MA = false    # toggle true unless you're sure VCF only has biallelic sites
        Boolean add_annotations = false    # toggle true if VCF missing SVTYPE, SVLEN, END fields for large INDELs
        Boolean perform_qc = true
        Boolean run_wittyer = true

        File? evaluation_bed
        Float? evaluation_pct    # Defaults to checking at least one base overlap evaluation_bed

        Array[File] bed_regions = []
        Array[String] bed_labels = []
        Int breakpoint_padding = 20

        # Types to stratify stats over for Truvari; possible values for SVTYPE in inputs INFO field
        # By default, workflow adds "ALL" category to run over entire VCF (subset first if only interested in one type)
        Array[String] sv_types = ["INS", "DEL", "DUP", "CNV", "INV", "CPX", "BND"]
    }

    Array[String] all_sv_types = flatten([["ALL"], sv_types])

    # Subset to evaluation regions
    if (defined(evaluation_bed)) {
        call SubsetEvaluation as SubsetComp {
            input:
                input_vcf=comp_vcf,
                input_vcf_index=comp_vcf_index,
                evaluation_bed=select_first([evaluation_bed]),
                evaluation_pct=evaluation_pct
        }

        if (defined(base_vcf)) {
            call SubsetEvaluation as SubsetBase {
                input:
                    input_vcf=select_first([base_vcf]),
                    input_vcf_index=select_first([base_vcf_index]),
                    evaluation_bed=select_first([evaluation_bed]),
                    evaluation_pct=evaluation_pct
            }
        }
    }

    if (split_MA) {
        call SplitMASites as SplitComp {
            input:
                input_vcf=select_first([SubsetComp.output_vcf, comp_vcf]),
                input_vcf_index=select_first([SubsetComp.output_vcf_index, comp_vcf_index])
        }

        if (defined(base_vcf)) {
            call SplitMASites as SplitBase {
                input:
                    input_vcf=select_first([SubsetBase.output_vcf, base_vcf]),
                    input_vcf_index=select_first([SubsetBase.output_vcf_index, base_vcf_index])
            }
        }
    }

    # Replace input files with subset versions
    File subset_comp_vcf = select_first([SplitComp.output_vcf, SubsetComp.output_vcf, comp_vcf])
    File subset_comp_vcf_index = select_first([SplitComp.output_vcf_index, SubsetComp.output_vcf_index, comp_vcf_index])
    File subset_base_vcf = select_first([SplitBase.output_vcf, SubsetBase.output_vcf, base_vcf])
    File subset_base_vcf_index = select_first([SplitBase.output_vcf_index, SubsetBase.output_vcf_index, base_vcf_index])


    # QC Tasks
    if (perform_qc) {
        call CollectQcMetrics as CompQcMetrics {
            input:
                input_vcf=subset_comp_vcf,
                input_vcf_index=subset_comp_vcf_index,
                experiment=experiment,
                label_suffix="Comp",
                ref_dict=ref_dict,
                bed_regions=bed_regions,
                bed_labels=bed_labels
        }

        if (defined(base_vcf)) {
            call CollectQcMetrics as BaseQcMetrics {
                input:
                    input_vcf=subset_base_vcf,
                    input_vcf_index=subset_base_vcf_index,
                    experiment=experiment,
                    label_suffix="Base",
                    ref_dict=ref_dict,
                    bed_regions=bed_regions,
                    bed_labels=bed_labels
            }

            call CombineSummaries as CombineQcMetrics {
                input:
                    tables=[CompQcMetrics.stats, BaseQcMetrics.stats],
                    output_file_name="combined_qc_stats.tsv"
            }
        }
    }


    # Benchmarking - Truvari tasks
    ## Run Truvari over different SV types for each sample name pair
    if (defined(base_vcf)) {
        scatter(sample_pair in zip(base_sample_names, comp_sample_names)) {
            scatter(sv_type in all_sv_types) {
                call RunTruvari {
                    input:
                        base_vcf=subset_base_vcf,
                        base_vcf_index=subset_base_vcf_index,
                        base_vcf_sample_name=sample_pair.left,
                        comp_vcf=subset_comp_vcf,
                        comp_vcf_index=subset_comp_vcf_index,
                        comp_vcf_sample_name=sample_pair.right,
                        experiment=experiment,
                        bed_regions=bed_regions,
                        bed_labels=bed_labels,
                        sv_type=sv_type,
                }
            }
            ## Collect Interval stats for Truvari FN/FP outputs
            call CollectIntervalComparisonMetrics {
                input:
                base_vcf=subset_base_vcf,
                base_vcf_index=subset_base_vcf_index,
                base_sample_name=sample_pair.left,
                comp_vcf=subset_comp_vcf,
                comp_vcf_index=subset_comp_vcf_index,
                comp_sample_name=sample_pair.right,
                tp_base=select_first(RunTruvari.tp_base),
                tp_base_index=select_first(RunTruvari.tp_base_index),
                tp_comp=select_first(RunTruvari.tp_comp),
                tp_comp_index=select_first(RunTruvari.tp_comp_index),
                fp=select_first(RunTruvari.fp),
                fp_index=select_first(RunTruvari.fp_index),
                fn=select_first(RunTruvari.fn),
                fn_index=select_first(RunTruvari.fn_index),
                bed_regions=bed_regions,
                bed_labels=bed_labels,
                breakpoint_padding=breakpoint_padding,
                ref_dict=ref_dict
            }
        }
    }

    ## Combine summaries across the shards
    if (defined(base_vcf)) {
        if (length(select_all(comp_sample_names)) > 0) {
            # First two: collect "bedtools closest" stats
            call CombineSummaries as CombineIntervalMetricsFPClosest {
                input:
                    tables=select_first([CollectIntervalComparisonMetrics.fp_closest]),
                    output_file_name="truvari_fp_closest.tsv"
            }

            call CombineSummaries as CombineIntervalMetricsFNClosest {
                input:
                    tables=select_first([CollectIntervalComparisonMetrics.fn_closest]),
                    output_file_name="truvari_fn_closest.tsv"
            }

            # Next four: collect interval overlap stats
            call CombineSummaries as CombineIntervalMetricsTPBaseIntervals {
                input:
                    tables=select_first([CollectIntervalComparisonMetrics.tp_base_intervals]),
                    output_file_name="truvari_tp-base_intervals.tsv"
            }

            call CombineSummaries as CombineIntervalMetricsTPCompIntervals {
                input:
                    tables=select_first([CollectIntervalComparisonMetrics.tp_comp_intervals]),
                    output_file_name="truvari_tp-comp_intervals.tsv"
            }

            call CombineSummaries as CombineIntervalMetricsFPIntervals {
                input:
                    tables=select_first([CollectIntervalComparisonMetrics.fp_intervals]),
                    output_file_name="truvari_fp_intervals.tsv"
            }

            call CombineSummaries as CombineIntervalMetricsFNIntervals {
                input:
                    tables=select_first([CollectIntervalComparisonMetrics.fn_intervals]),
                    output_file_name="truvari_fn_intervals.tsv"
            }

            # Finally: collect all basic Truvari stats
            call CombineSummaries as CombineBenchStats {
                input:
                    tables=select_all(flatten(select_first([RunTruvari.summary]))),
                    output_file_name="truvari_bench_summary.tsv"
            }
        }
    }


    # Benchmarking - Wittyer tasks
    if (defined(base_vcf) && run_wittyer) {
        if (length(select_all(comp_sample_names)) > 0) {
            call WittyerEval as Wittyer {
                input:
                    truth_vcf=subset_base_vcf,
                    eval_vcf=subset_comp_vcf,
                    truth_sample_names=base_sample_names,
                    query_sample_names=comp_sample_names
            }

            call CleanBasicWittyerStats {
                input:
                    wittyer_stats=Wittyer.wittyer_stats,
                    experiment=experiment
            }

            scatter(sample_pair in zip(Wittyer.wittyer_annotated_vcf, zip(select_first([base_sample_names]), select_first([comp_sample_names])))) {
                call AddIntervalOverlapStatsWittyer as WittyerIntervals {
                    input:
                        input_vcf=sample_pair.left,
                        base_sample_name=sample_pair.right.left,
                        comp_sample_name=sample_pair.right.right,
                        experiment=experiment,
                        ref_dict=ref_dict,
                        bed_regions=select_all(bed_regions),
                        bed_labels=select_all(bed_labels),
                        breakpoint_padding=breakpoint_padding
                }
            }

            call CombineSummaries as CombineWittierTruthOverlapStats {
                input:
                    tables=WittyerIntervals.truth_table_intervals,
                    output_file_name="wittyer_truth_intervals.tsv"
            }

            call CombineSummaries as CombineWittierQueryOverlapStats {
                input:
                    tables=WittyerIntervals.query_table_intervals,
                    output_file_name="wittyer_query_intervals.tsv"
            }

            call CombineSummaries as CombineWittierNoGtOverlapStats {
                input:
                    tables=WittyerIntervals.no_gt_intervals,
                    output_file_name="wittyer_nogt_intervals.tsv"
            }
        }
    }

    call CombineFiles {
        input:
            qc_files=[select_first([CombineQcMetrics.combined_file, CompQcMetrics.stats])],
            truvari_files=select_all([CombineBenchStats.combined_file, CombineIntervalMetricsFPClosest.combined_file, CombineIntervalMetricsFNClosest.combined_file,
                CombineIntervalMetricsTPBaseIntervals.combined_file, CombineIntervalMetricsTPCompIntervals.combined_file, CombineIntervalMetricsFPIntervals.combined_file,
                CombineIntervalMetricsFNIntervals.combined_file]),
            wittyer_files=select_all([CleanBasicWittyerStats.wittyer_stats_cleaned, CombineWittierTruthOverlapStats.combined_file, CombineWittierQueryOverlapStats.combined_file,
                CombineWittierNoGtOverlapStats.combined_file])
    }

    output {
        # QC outputs
        File? qc_summary = select_first([CombineQcMetrics.combined_file, CompQcMetrics.stats])

        # Truvari outputs
        File? truvari_bench_summary = CombineBenchStats.combined_file
        File? truvari_fp_closest = CombineIntervalMetricsFPClosest.combined_file
        File? truvari_fn_closest = CombineIntervalMetricsFNClosest.combined_file
        File? truvari_tp_base_intervals = CombineIntervalMetricsTPBaseIntervals.combined_file
        File? truvari_tp_comp_intervals = CombineIntervalMetricsTPCompIntervals.combined_file
        File? truvari_fp_intervals = CombineIntervalMetricsFPIntervals.combined_file
        File? truvari_fn_intervals = CombineIntervalMetricsFNIntervals.combined_file

        # Wittyer outputs
        File? wittyer_basic_stats = CleanBasicWittyerStats.wittyer_stats_cleaned
        File? wittyer_truth_intervals = CombineWittierTruthOverlapStats.combined_file
        File? wittyer_query_intervals = CombineWittierQueryOverlapStats.combined_file
        File? wittyer_nogt_intervals = CombineWittierNoGtOverlapStats.combined_file

        # Compressed/combined outputs
        File? combined_files = CombineFiles.combined_files
    }
}

task SubsetEvaluation {
    input {
        File input_vcf
        File input_vcf_index

        File evaluation_bed
        Float? evaluation_pct

        # Runtime parameters
        Int disk_size = ceil(2 * size(input_vcf, "GB")) + 100
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xueo pipefail

        # Subset input VCF to only sites which intersect evaluation regions over pct overlap threshold
        bcftools view -h ~{input_vcf} > header.txt
        bedtools intersect -a ~{input_vcf} -b ~{evaluation_bed} ~{"-f" + evaluation_pct} -u > variants.vcf

        cat header.txt variants.vcf | gzip > output.vcf.gz
        bcftools index -t output.vcf.gz

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.3"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
    }

    output {
        File output_vcf = "output.vcf.gz"
        File output_vcf_index = "output.vcf.gz.tbi"
    }
}

task SplitMASites {
    input {
        File input_vcf
        File input_vcf_index

        # Runtime parameters
        Int disk_size = ceil(2 * size(input_vcf, "GB")) + 100
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xueo pipefail

        bcftools norm -m -any ~{input_vcf} > split.vcf.gz
        bcftools index -t split.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.3"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
    }

    output {
        File output_vcf = "split.vcf.gz"
        File output_vcf_index = "split.vcf.gz.tbi"
    }
}

task CollectQcMetrics {
    input {
        File input_vcf
        File input_vcf_index

        String experiment = ""
        String label_suffix = "Comp"

        File? ref_dict

        Array[File] bed_regions = []
        Array[String] bed_labels = []
        Int? breakpoint_padding = 20

        # Runtime parameters
        Int disk_size = ceil(2 * size(input_vcf, "GB")) + 100
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xueo pipefail

        # Add allele count annotations, in case they're missing
        bcftools +fill-tags -o tagged.vcf.gz "~{input_vcf}"
        bcftools index -t -f tagged.vcf.gz

        # Main query formatting
        MAIN_HEADER="CHROM\tPOS\tEND\tFILTER\tSVTYPE\tSVLEN\tAC\tAC_Het\tAC_Hom\tAN\tExcHet\tHWE\tAF\tMAF\tNS\tExperiment"
        MAIN_QUERY="%CHROM\t%POS\t%INFO/END\t%FILTER\t%INFO/SVTYPE\t%INFO/SVLEN\t%AC\t%AC_Het\t%AC_Hom\t%AN\t%ExcHet\t%HWE\t%AF\t%MAF\t%NS"

        # Add newline to query expression
        MAIN_QUERY="${MAIN_QUERY}\n"

        # Create table of relevant stats from variant annotations
        # Requires input to have INFO fields: END, SVTYPE, SVLEN
        echo -e "${MAIN_HEADER}" > qc_stats.tsv
        bcftools query \
            -f"${MAIN_QUERY}" \
            tagged.vcf.gz | awk -v OFS='\t' '{ print $0, "~{experiment}-~{label_suffix}" }' >> qc_stats-pre_intervals.tsv

        # Add interval overlap data
        # Clean ref_dict into genome file for bedtools
        tail -n +2 ~{ref_dict} | cut -f2,3 - | sed 's/SN://g' - | sed 's/LN://g' - > ref.genome

        # Collect interval bed overlap stats
        generate_interval_stats () {
            CURRENT_LABEL=$1
            INTERVAL_FILE=$2
            INPUT_LABEL=$3
            VCF_FILE=$4

            echo -e "${CURRENT_LABEL}-count\t${CURRENT_LABEL}-overlap" > header.txt
            # WARNING: annotate does not preserve order of input bed regions!
            cut -f1-3 $VCF_FILE > vcf.bed
            bedtools annotate -both -i vcf.bed -files $INTERVAL_FILE | bedtools sort -i - | rev | cut -f1-2 | rev > "${CURRENT_LABEL}-${INPUT_LABEL}-preheader.bed"
            cat header.txt "${CURRENT_LABEL}-${INPUT_LABEL}-preheader.bed" > "${CURRENT_LABEL}-${INPUT_LABEL}-annotated.bed"

            # Also add in breakpoint overlap stats using POS and END values
            echo -e "${CURRENT_LABEL}-LBEND-count\t${CURRENT_LABEL}-LBEND-overlap" > pos-header.txt
            echo -e "${CURRENT_LABEL}-RBEND-count\t${CURRENT_LABEL}-RBEND-overlap" > end-header.txt
            awk '{OFS="\t"; print $1, $2, $2}' vcf.bed | bedtools slop -b ~{breakpoint_padding} -i - -g ref.genome > pos.bed
            awk '{OFS="\t"; print $1, $3, $3}' vcf.bed | bedtools slop -b ~{breakpoint_padding} -i - -g ref.genome > end.bed
            bedtools annotate -both -i pos.bed -files $INTERVAL_FILE | bedtools sort -i - | rev | cut -f1-2 | rev > "${CURRENT_LABEL}-${INPUT_LABEL}-pos-preheader.bed"
            bedtools annotate -both -i end.bed -files $INTERVAL_FILE | bedtools sort -i - | rev | cut -f1-2 | rev > "${CURRENT_LABEL}-${INPUT_LABEL}-end-preheader.bed"
            cat pos-header.txt "${CURRENT_LABEL}-${INPUT_LABEL}-pos-preheader.bed" > "${CURRENT_LABEL}-${INPUT_LABEL}-pos-annotated.bed"
            cat end-header.txt "${CURRENT_LABEL}-${INPUT_LABEL}-end-preheader.bed" > "${CURRENT_LABEL}-${INPUT_LABEL}-end-annotated.bed"

            # Paste together three overlap stats files
            paste "${CURRENT_LABEL}-${INPUT_LABEL}-annotated.bed" "${CURRENT_LABEL}-${INPUT_LABEL}-pos-annotated.bed" "${CURRENT_LABEL}-${INPUT_LABEL}-end-annotated.bed" > "${CURRENT_LABEL}-${INPUT_LABEL}-full-annotated.bed"
        }

        if [ ~{length(bed_regions)} -gt 0 ]
        then
            INTERVAL_FILES=(~{sep=' ' bed_regions})
            INTERVAL_LABELS=(~{sep=' ' bed_labels})

            for i in {0..~{length(bed_regions)-1}}
            do
                CURRENT_LABEL="${INTERVAL_LABELS[$i]}"
                CURRENT_FILE="${INTERVAL_FILES[$i]}"
                generate_interval_stats $CURRENT_LABEL $CURRENT_FILE "qc" qc_stats-pre_intervals.tsv
            done

            # Combine across interval beds
            echo -e $MAIN_HEADER > header.txt
            cat header.txt qc_stats-pre_intervals.tsv > qc_stats-header.tsv
            paste qc_stats-header.tsv *-qc-full-annotated.bed > qc_stats.tsv
        else
            mv qc_stats-pre_header.tsv qc_stats.tsv
        fi

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.3"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
    }

    output {
        File stats = "qc_stats.tsv"
    }
}

task RunTruvari {
    input {
        File base_vcf
        File base_vcf_index
        String? base_vcf_sample_name

        File comp_vcf
        File comp_vcf_index
        String? comp_vcf_sample_name

        String experiment = "Experiment"
        String sv_type

        # Tool arguments
        ## Comparison Threshold
        Int ref_dist = 500
        Float pct_seq = 0.7
        Int min_hap_len = 50
        Float pct_size = 0.7
        Float pct_overlap = 0.0
        Boolean type_ignore = false    # Toggle true to ignore matching SV types for TP
        Boolean dup_to_ins = false    # Toggle true to treat DUP as INS
        Int chunk_size = 1000    # Max ref distance to compare calls
        String num_matches = "single"   # Num times a variant can match; single, ac, or multi

        ## Filtering
        Int call_size_min = 50
        Int base_size_min = 50
        Int size_max = 50000
        Boolean pass_only = true
        Array[File] bed_regions = []    # Restrict to comparing SVs overlapping input beds
        Array[String] bed_labels = []
        Int extend = 0     # Padding to add onto bed_region

        Boolean debug_mode = false

        # Runtime parameters
        Int disk_size = ceil(2 * size(base_vcf, "GB") + 2 * size(comp_vcf, "GB")) + 100
        Int cpu = 8
        Int memory_ram = 16
    }

    command <<<
        set -xueo pipefail

        # Subset VCFs by SV type
        if [ "~{sv_type}" = "ALL" ]
        then
            ln -s ~{base_vcf} base-subset.vcf.gz
            ln -s ~{base_vcf_index} base-subset.vcf.gz.tbi
            ln -s ~{comp_vcf} comp-subset.vcf.gz
            ln -s ~{comp_vcf_index} comp-subset.vcf.gz.tbi
        else
            bcftools view -i 'INFO/SVTYPE == "~{sv_type}"' -o base-subset.vcf.gz ~{base_vcf}
            bcftools index -t -f base-subset.vcf.gz
            bcftools view -i 'INFO/SVTYPE == "~{sv_type}"' -o comp-subset.vcf.gz ~{comp_vcf}
            bcftools index -t -f comp-subset.vcf.gz
        fi

        # Subset to sites for given sample names since Truvari flags don't seem to work right
        if [ -z "~{base_vcf_sample_name}~{comp_vcf_sample_name}" ]
        then
            mv base-subset.vcf.gz base.vcf.gz
            mv base-subset.vcf.gz.tbi base.vcf.gz.tbi
            mv comp-subset.vcf.gz comp.vcf.gz
            mv comp-subset.vcf.gz.tbi comp.vcf.gz.tbi
        else
            bcftools view ~{"-s" + base_vcf_sample_name} --min-ac 1 -o base.vcf.gz base-subset.vcf.gz
            bcftools index -t base.vcf.gz
            bcftools view ~{"-s" + comp_vcf_sample_name} --min-ac 1 -o comp.vcf.gz comp-subset.vcf.gz
            bcftools index -t comp.vcf.gz
        fi


        # Function for running Truvari over a subset cut out from intervals
        run_truvari_on_subset() {
            # Parse input variables about intervals
            INTERVAL_FILE=$1
            INTERVAL_LABEL=$2

            # Determine whether to use includebed argument based on whether interval file is provided
            INCLUDE_BED_COMMAND=""
            if [ "$INTERVAL_FILE" = " " ]
            then
                INCLUDE_BED_COMMAND=""
            else
                INCLUDE_BED_COMMAND="--includebed ${INTERVAL_FILE}"
            fi

            # Run Truvari for benchmarking
            truvari bench \
                -b base.vcf.gz \
                ~{"--bSample " + base_vcf_sample_name} \
                -c comp.vcf.gz \
                ~{"--cSample " + comp_vcf_sample_name} \
                -o "output_dir-${INTERVAL_LABEL}" \
                ~{true="--debug" false="" debug_mode} \
                $INCLUDE_BED_COMMAND \
                ~{"--refdist " + ref_dist} \
                ~{"--pctseq " + pct_seq} \
                ~{"--minhaplen " + min_hap_len} \
                ~{"--pctsize " + pct_size} \
                ~{"--pctovl " + pct_overlap} \
                ~{true="--typeignore " false="" type_ignore} \
                ~{true="--dup-to-ins " false="" dup_to_ins} \
                ~{"--chunksize " + chunk_size} \
                ~{"--sizemin " + call_size_min} \
                ~{"--sizefilt " + base_size_min} \
                ~{"--sizemax " + size_max} \
                ~{true="--passonly " false="" pass_only} \
                ~{"--extend " + extend} \
                ~{"--pick " + num_matches}
        }

        # Run Truvari over all interval files, SV types, and evidence classes
        INTERVAL_FILES=(" " ~{sep=' ' bed_regions})
        INTERVAL_LABELS=("WholeGenome" ~{sep=' ' bed_labels})

        for i in {0..~{length(bed_regions)}}
        do
            run_truvari_on_subset "${INTERVAL_FILES[$i]}" "${INTERVAL_LABELS[$i]}"
        done

        # Use Python to collect the output results and compile across intervals into one table
        python3 << CODE

        import pandas as pd
        import json

        # Combine Truvari outputs on TP, FP, etc. across different subsets into table
        full_df = pd.DataFrame()
        for interval in ["WholeGenome"] + ["~{default="" sep="\", \"" bed_labels}"]:
            with open(f'output_dir-{interval}/summary.json') as file:
                json_file = json.load(file)

            df = pd.DataFrame({k: v for k,v in json_file.items() if k != 'gt_matrix'}, index=[0])
            df['Base_Name'] = "~{default="" base_vcf_sample_name}"
            df['Comp_Name'] = "~{default="" comp_vcf_sample_name}"
            df['Experiment'] = "~{experiment}"
            df['Interval'] = interval
            df['SV_Type'] = "~{sv_type}"

            full_df = pd.concat([full_df, df])

        full_df.to_csv("CombinedSummary.tsv", sep='\t', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.3"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
    }

    output {
        File fn = "output_dir-WholeGenome/fn.vcf.gz"
        File fn_index = "output_dir-WholeGenome/fn.vcf.gz.tbi"
        File fp = "output_dir-WholeGenome/fp.vcf.gz"
        File fp_index = "output_dir-WholeGenome/fp.vcf.gz.tbi"
        File tp_base = "output_dir-WholeGenome/tp-base.vcf.gz"
        File tp_base_index = "output_dir-WholeGenome/tp-base.vcf.gz.tbi"
        File tp_comp = "output_dir-WholeGenome/tp-comp.vcf.gz"
        File tp_comp_index = "output_dir-WholeGenome/tp-comp.vcf.gz.tbi"

        File summary = "CombinedSummary.tsv"
    }
}

task CollectIntervalComparisonMetrics {
    input {
        File base_vcf
        File base_vcf_index
        String? base_sample_name

        File comp_vcf
        File comp_vcf_index
        String? comp_sample_name

        File tp_base
        File tp_base_index
        File tp_comp
        File tp_comp_index
        File fp
        File fp_index
        File fn
        File fn_index

        Array[File] bed_regions = []
        Array[String] bed_labels = []
        Int breakpoint_padding = 20
        File ref_dict

        Int k_closest = 3    # Number of close by variants to compare to

        # Runtime parameters
        Int disk_size = ceil(2 * size(comp_vcf, "GB")) + 100
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xueo pipefail

        # Transform VCFs into bed files for more convenient analysis
        bcftools view ~{"-s" + base_sample_name} --min-ac 1 "~{base_vcf}" | bcftools query -f'%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t%FILTER\n' - | bedtools sort -i - > base.bed
        bcftools view ~{"-s" + base_sample_name} --min-ac 1 "~{tp_base}" | bcftools query -f'%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t%FILTER\n' - | bedtools sort -i - > tp-base.bed
        bcftools view ~{"-s" + base_sample_name} --min-ac 1 "~{fn}" | bcftools query -f'%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t%FILTER\n' - | bedtools sort -i - > fn.bed
        bcftools view ~{"-s" + comp_sample_name} --min-ac 1 "~{comp_vcf}" | bcftools query -f'%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t%FILTER\n' - | bedtools sort -i - > comp.bed
        bcftools view ~{"-s" + comp_sample_name} --min-ac 1 "~{tp_comp}" | bcftools query -f'%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t%FILTER\n' - | bedtools sort -i - > tp-comp.bed
        bcftools view ~{"-s" + comp_sample_name} --min-ac 1 "~{fp}" | bcftools query -f'%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t%FILTER\n' - | bedtools sort -i - > fp.bed

        # Generate "closest" data
        echo -e "BASENAME\tCOMPNAME\tLCHROM\tLPOS\tLEND\tLLEN\tLTYPE\tLFILTER\tRCHROM\tRPOS\tREND\tRLEN\tRTYPE\tRFILTER\tDIST" > header.txt

        bedtools closest -a fp.bed -b base.bed -k ~{k_closest} -D ref > fp-closest.bed
        awk -v OFS='\t' '{ print "~{base_sample_name}", "~{comp_sample_name}", $0 }' fp-closest.bed > fp-closest-samples.bed
        cat header.txt fp-closest-samples.bed > fp-closest-final.bed

        bedtools closest -a fn.bed -b comp.bed -k ~{k_closest} -D ref > fn-closest.bed
        awk -v OFS='\t' '{ print "~{base_sample_name}", "~{comp_sample_name}", $0 }' fn-closest.bed > fn-closest-samples.bed
        cat header.txt fn-closest-samples.bed > fn-closest-final.bed

        # Clean ref_dict into genome file for bedtools
        tail -n +2 ~{ref_dict} | cut -f2,3 - | sed 's/SN://g' - | sed 's/LN://g' - > ref.genome

        # Collect interval bed overlap stats
        generate_interval_stats () {
            CURRENT_LABEL=$1
            INTERVAL_FILE=$2
            INPUT_LABEL=$3
            VCF_FILE=$4

            echo -e "${CURRENT_LABEL}-count\t${CURRENT_LABEL}-overlap" > header.txt
            # WARNING: annotate does not preserve order of input bed regions!
            cut -f1-3 $VCF_FILE > vcf.bed
            bedtools annotate -both -i vcf.bed -files $INTERVAL_FILE | bedtools sort -i - | rev | cut -f1-2 | rev > "${CURRENT_LABEL}-${INPUT_LABEL}-preheader.bed"
            cat header.txt "${CURRENT_LABEL}-${INPUT_LABEL}-preheader.bed" > "${CURRENT_LABEL}-${INPUT_LABEL}-annotated.bed"

            # Also add in breakpoint overlap stats using POS and END values
            echo -e "${CURRENT_LABEL}-LBEND-count\t${CURRENT_LABEL}-LBEND-overlap" > pos-header.txt
            echo -e "${CURRENT_LABEL}-RBEND-count\t${CURRENT_LABEL}-RBEND-overlap" > end-header.txt
            awk '{OFS="\t"; print $1, $2, $2}' vcf.bed | bedtools slop -b ~{breakpoint_padding} -i - -g ref.genome > pos.bed
            awk '{OFS="\t"; print $1, $3, $3}' vcf.bed | bedtools slop -b ~{breakpoint_padding} -i - -g ref.genome > end.bed
            bedtools annotate -both -i pos.bed -files $INTERVAL_FILE | bedtools sort -i - | rev | cut -f1-2 | rev > "${CURRENT_LABEL}-${INPUT_LABEL}-pos-preheader.bed"
            bedtools annotate -both -i end.bed -files $INTERVAL_FILE | bedtools sort -i - | rev | cut -f1-2 | rev > "${CURRENT_LABEL}-${INPUT_LABEL}-end-preheader.bed"
            cat pos-header.txt "${CURRENT_LABEL}-${INPUT_LABEL}-pos-preheader.bed" > "${CURRENT_LABEL}-${INPUT_LABEL}-pos-annotated.bed"
            cat end-header.txt "${CURRENT_LABEL}-${INPUT_LABEL}-end-preheader.bed" > "${CURRENT_LABEL}-${INPUT_LABEL}-end-annotated.bed"

            # Paste together three overlap stats files
            paste "${CURRENT_LABEL}-${INPUT_LABEL}-annotated.bed" "${CURRENT_LABEL}-${INPUT_LABEL}-pos-annotated.bed" "${CURRENT_LABEL}-${INPUT_LABEL}-end-annotated.bed" > "${CURRENT_LABEL}-${INPUT_LABEL}-full-annotated.bed"
        }

        if [ ~{length(bed_regions)} -gt 0 ]
        then
            INTERVAL_FILES=(~{sep=' ' bed_regions})
            INTERVAL_LABELS=(~{sep=' ' bed_labels})

            for i in {0..~{length(bed_regions)-1}}
            do
                CURRENT_LABEL="${INTERVAL_LABELS[$i]}"
                CURRENT_FILE="${INTERVAL_FILES[$i]}"
                generate_interval_stats $CURRENT_LABEL $CURRENT_FILE "tp-base" tp-base.bed
                generate_interval_stats $CURRENT_LABEL $CURRENT_FILE "tp-comp" tp-comp.bed
                generate_interval_stats $CURRENT_LABEL $CURRENT_FILE "fp" fp.bed
                generate_interval_stats $CURRENT_LABEL $CURRENT_FILE "fn" fn.bed
            done

            # Combine across interval beds and add sample names
            echo -e "BASENAME\tCOMPNAME\tCHROM\tPOS\tEND\tSVLEN\tSVTYPE\tFILTER" > header.txt

            awk -v OFS='\t' '{ print "~{base_sample_name}", "~{comp_sample_name}", $0 }' tp-base.bed > tp-base_samples.bed
            cat header.txt tp-base_samples.bed > tp-base_header.bed
            paste tp-base_header.bed *-tp-base-full-annotated.bed > tp-base_intervals-final.bed

            awk -v OFS='\t' '{ print "~{base_sample_name}", "~{comp_sample_name}", $0 }' tp-comp.bed > tp-comp_samples.bed
            cat header.txt tp-comp_samples.bed > tp-comp_header.bed
            paste tp-comp_header.bed *-tp-comp-full-annotated.bed > tp-comp_intervals-final.bed

            awk -v OFS='\t' '{ print "~{base_sample_name}", "~{comp_sample_name}", $0 }' fp.bed > fp_samples.bed
            cat header.txt fp_samples.bed > fp_header.bed
            paste fp_header.bed *-fp-full-annotated.bed > fp_intervals-final.bed

            awk -v OFS='\t' '{ print "~{base_sample_name}", "~{comp_sample_name}", $0 }' fn.bed > fn_samples.bed
            cat header.txt fn_samples.bed > fn_header.bed
            paste fn_header.bed *-fn-full-annotated.bed > fn_intervals-final.bed
        fi

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.3"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
    }

    output {
        File fp_closest = "fp-closest-final.bed"
        File fn_closest = "fn-closest-final.bed"

        File tp_base_intervals = "tp-base_intervals-final.bed"
        File tp_comp_intervals = "tp-comp_intervals-final.bed"
        File fp_intervals = "fp_intervals-final.bed"
        File fn_intervals = "fn_intervals-final.bed"
    }
}

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
        String wittyer_docker = "us.gcr.io/broad-dsde-methods/wittyer:v1.0"

        Int? mem
        Int? disk_space
    }

    # If mem and disk size were not specified, use 4GB and 100 GB as default
    Int mem_size = select_first([mem, 4])
    Int disk_size = select_first([disk_space,100])

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
        docker: wittyer_docker
        bootDiskSizeGb: 12
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
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

        File ref_dict

        # Runtime parameters
        Int disk_size = ceil(2 * size(input_vcf, "GB")) + 100
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xueo pipefail

        # Extract useful columns from Wittyer VCF
        bcftools query -f'%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t%FILTER\t%WHERE\t%WIN[\t%GT\t%WIT\t%WHY]\n' ~{input_vcf} | bedtools sort -i - > wittyer_vcf_labels.tsv

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

        # Clean ref_dict into genome file for bedtools
        tail -n +2 ~{ref_dict} | cut -f2,3 - | sed 's/SN://g' - | sed 's/LN://g' - > ref.genome

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

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.3"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
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
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.3"
        disks: "local-disk " + 50 + " HDD"
        memory: 4 + " GB"
        cpu: 2
    }

    output {
        File wittyer_stats_cleaned = "wittyer_stats-cleaned.tsv"
    }
}

task CombineSummaries {
    input {
        Array[File] tables
        String output_file_name

        # Runtime parameters
        Int disk_size = ceil(2 * size(tables, "GB")) + 50
        Int cpu = 4
        Int memory_ram = 8
    }

    command <<<
        set -xueo pipefail

        python3 << CODE

        import pandas as pd

        full_df = pd.DataFrame()
        for file in ["~{default="" sep="\", \"" tables}"]:
            df = pd.read_csv(file, sep='\t')
            full_df = pd.concat([full_df, df])

        full_df.to_csv('~{output_file_name}', sep='\t', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.3"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
    }

    output {
        File combined_file = "~{output_file_name}"
    }
}

task CombineFiles {
    input {
        Array[File] qc_files
        Array[File] truvari_files
        Array[File] wittyer_files
    }

    command <<<
        set -xueo pipefail

        mkdir benchmark_outputs

        mkdir benchmark_outputs/qc_files
        cp ~{sep=" " qc_files} benchmark_outputs/qc_files/

        mkdir benchmark_outputs/truvari_files
        cp ~{sep=" " truvari_files} benchmark_outputs/truvari_files/

        mkdir benchmark_outputs/wittyer_files
        cp ~{sep=" " wittyer_files} benchmark_outputs/wittyer_files/

        tar -zcvf benchmark_outputs.tar.gz benchmark_outputs
    >>>

    runtime {
        docker: "ubuntu:23.10"
        disks: "local-disk 250 HDD"
        memory: "16GB"
        cpu: 2
    }

    output {
        File combined_files = "benchmark_outputs.tar.gz"
    }
}