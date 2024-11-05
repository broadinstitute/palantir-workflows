version 1.0

import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/Utilities/WDLs/CreateIGVSession.wdl" as IGV
import "tasks.wdl" as SV_tasks

struct RuntimeAttributes {
    Int disk_size
    Int cpu
    Int memory
}

workflow BenchmarkSVs {
    input {
        File base_vcf
        File base_vcf_index
        String base_sample_name

        File comp_vcf
        File comp_vcf_index
        String comp_sample_name

        String? experiment

        File ref_fasta
        File ref_fai

        File? evaluation_bed
        Float? evaluation_pct    # Defaults to checking at least one base overlap evaluation_bed; between 0 and 1

        Array[File] bed_regions = []
        Array[String] bed_labels = []
        Int breakpoint_padding = 20

        Array[Int] svlen_bin_cutoffs = [100, 250, 1000, 2500, 10000, 25000]

        Boolean create_igv_session = true
        Array[File]? optional_igv_bams
    }

    # Subset to evaluation regions for remainder of analysis
    if (defined(evaluation_bed)) {
        call SV_tasks.SubsetEvaluation as SubsetComp {
            input:
                input_vcf=comp_vcf,
                input_vcf_index=comp_vcf_index,
                evaluation_bed=select_first([evaluation_bed]),
                evaluation_pct=evaluation_pct
        }

        if (defined(base_vcf)) {
            call SV_tasks.SubsetEvaluation as SubsetBase {
                input:
                    input_vcf=select_first([base_vcf]),
                    input_vcf_index=select_first([base_vcf_index]),
                    evaluation_bed=select_first([evaluation_bed]),
                    evaluation_pct=evaluation_pct
            }
        }
    }

    # Replace input files with subset versions
    File subset_comp_vcf = select_first([SubsetComp.output_vcf, comp_vcf])
    File subset_comp_vcf_index = select_first([SubsetComp.output_vcf_index, comp_vcf_index])
    File subset_base_vcf = select_first([SubsetBase.output_vcf, base_vcf])
    File subset_base_vcf_index = select_first([SubsetBase.output_vcf_index, base_vcf_index])

    # Benchmarking - Truvari tasks
    ## Run Truvari over each sample name pair
    call RunTruvari {
        input:
            base_vcf=subset_base_vcf,
            base_vcf_index=subset_base_vcf_index,
            base_sample_name=base_sample_name,
            comp_vcf=subset_comp_vcf,
            comp_vcf_index=subset_comp_vcf_index,
            comp_sample_name=comp_sample_name,
            experiment=experiment,
            ref_fasta=ref_fasta,
            ref_fai=ref_fai
    }

    call SV_tasks.AddIntervalOverlapStats as TruvariTpBaseIntervals {
        input:
            input_table=RunTruvari.tp_base_table,
            table_label="truvari_bench",
            output_name=base_sample_name+"-tp_base_intervals",
            bed_regions=bed_regions,
            bed_labels=bed_labels,
            ref_fai=ref_fai,
            breakpoint_padding=breakpoint_padding
    }

    call SV_tasks.AddIntervalOverlapStats as TruvariTpCompIntervals {
        input:
            input_table=RunTruvari.tp_comp_table,
            table_label="truvari_bench",
            output_name=comp_sample_name+"-tp_comp_intervals",
            bed_regions=bed_regions,
            bed_labels=bed_labels,
            ref_fai=ref_fai,
            breakpoint_padding=breakpoint_padding
    }

    call SV_tasks.AddIntervalOverlapStats as TruvariFpIntervals {
        input:
            input_table=RunTruvari.fp_table,
            table_label="truvari_bench",
            output_name=comp_sample_name+"-fp_intervals",
            bed_regions=bed_regions,
            bed_labels=bed_labels,
            ref_fai=ref_fai,
            breakpoint_padding=breakpoint_padding
    }

    call SV_tasks.AddIntervalOverlapStats as TruvariFnIntervals {
        input:
            input_table=RunTruvari.fn_table,
            table_label="truvari_bench",
            output_name=base_sample_name+"-fn_intervals",
            bed_regions=bed_regions,
            bed_labels=bed_labels,
            ref_fai=ref_fai,
            breakpoint_padding=breakpoint_padding
    }

    call ComputeTruvariIntervalSummaryStats {
        input:
            base_sample_name=base_sample_name,
            comp_sample_name=comp_sample_name,
            tp_base_intervals=TruvariTpBaseIntervals.output_table,
            tp_comp_intervals=TruvariTpCompIntervals.output_table,
            fp_intervals=TruvariFpIntervals.output_table,
            fn_intervals=TruvariFnIntervals.output_table,
            bed_labels=bed_labels,
            experiment=experiment,
            svlen_bin_cutoffs=svlen_bin_cutoffs
    }

    call CollectTruvariClosestStats {
        input:
            base_table=RunTruvari.base_table,
            comp_table=RunTruvari.comp_table,
            fp_table=RunTruvari.fp_table,
            fn_table=RunTruvari.fn_table,
            base_sample_name=base_sample_name,
            comp_sample_name=comp_sample_name,
            experiment=experiment,
            ref_fai=ref_fai
    }

    if (create_igv_session) {
        call IGV.CreateIGVSession as IGVSession {
            input:
                bams=optional_igv_bams,
                vcfs=[RunTruvari.tp_base, RunTruvari.tp_comp, RunTruvari.fp, RunTruvari.fn],
                interval_files=bed_regions,
                reference=ref_fasta
        }
    }

    output {
        # Truvari outputs
        File truvari_bench_summary = ComputeTruvariIntervalSummaryStats.full_summary
        File? truvari_fp_closest = CollectTruvariClosestStats.fp_closest
        File? truvari_fn_closest = CollectTruvariClosestStats.fn_closest

        # IGV Session
         File? igv_session = IGVSession.igv_session
    }
}

task RunTruvari {
    input {
        File base_vcf
        File base_vcf_index
        String? base_sample_name

        File comp_vcf
        File comp_vcf_index
        String? comp_sample_name

        String experiment = "Experiment"

        File ref_fasta
        File ref_fai

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
        String num_matches = "ac"   # Num times a variant can match; single, ac, or multi

        ## Filtering
        Int call_size_min = 50
        Int base_size_min = 50
        Int size_max = 50000
        Boolean pass_only = true
        Array[String] drop_fields = ["FMT/AD"]    # For mismatching headers in refine step

        ## Refine args
        Boolean harmonize_phased_variants = false    # Use truvari refine on phased inputs
        String aligner = "mafft"

        Boolean debug_mode = false

        # Runtime parameters
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(base_vcf, "GB") + 2 * size(comp_vcf, "GB")) + 100, "cpu": 8, "memory": 16}
    }

    command <<<
        set -xueo pipefail

        # Subset to sites for given sample names since Truvari flags don't seem to work right
        if [ -z "~{base_sample_name}~{comp_sample_name}" ]
        then
            mv ~{base_vcf} base-mac1.vcf.gz
            mv ~{base_vcf_index} base-mac1.vcf.gz.tbi
            mv ~{comp_vcf} comp-mac1.vcf.gz
            mv ~{comp_vcf_index} comp-mac1.vcf.gz.tbi
        else
            bcftools view ~{"-s " + base_sample_name} --min-ac 1 -o base-mac1.vcf.gz -Wtbi ~{base_vcf}
            bcftools view ~{"-s " + comp_sample_name} --min-ac 1 -o comp-mac1.vcf.gz -Wtbi ~{comp_vcf}
        fi

        bcftools annotate -x ~{sep="," drop_fields} -o base.vcf.gz -Wtbi base-mac1.vcf.gz
        bcftools annotate -x ~{sep="," drop_fields} -o comp.vcf.gz -Wtbi comp-mac1.vcf.gz

        # Run Truvari for benchmarking
        truvari bench \
            -b base.vcf.gz \
            ~{"--bSample " + base_sample_name} \
            -c comp.vcf.gz \
            ~{"--cSample " + comp_sample_name} \
            -o "output_dir" \
            ~{true="--debug" false="" debug_mode} \
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
            ~{"--pick " + num_matches}

        RESULTS_STEM="output_dir/"
        if [ ~{harmonize_phased_variants} = true ]; then
            truvari refine \
                --reference ~{ref_fasta} \
                --regions output_dir/candidate.refine.bed \
                --recount \
                --use-region-coords \
                --use-original-vcfs \
                --align ~{aligner} \
                outputs

            truvari anno svinfo output_dir/phab_bench/tp-base.vcf.gz -o output_dir/phab_bench/anno-tp-base.vcf.gz
            bcftools index -t output_dir/phab_bench/anno-tp-base.vcf.gz
            truvari anno svinfo output_dir/phab_bench/tp-comp.vcf.gz -o output_dir/phab_bench/anno-tp-comp.vcf.gz
            bcftools index -t output_dir/phab_bench/anno-tp-comp.vcf.gz
            truvari anno svinfo output_dir/phab_bench/fp.vcf.gz -o output_dir/phab_bench/anno-fp.vcf.gz
            bcftools index -t output_dir/phab_bench/anno-fp.vcf.gz
            truvari anno svinfo output_dir/phab_bench/fn.vcf.gz -o output_dir/phab_bench/anno-fn.vcf.gz
            bcftools index -t output_dir/phab_bench/anno-fn.vcf.gz

            RESULTS_STEM="output_dir/phab_bench/anno-"
        fi

        # Use Python to collect the output results and compile across intervals into one table
        python3 << CODE
        import pandas as pd
        import json

        output_dir = 'output_dir' if "~{harmonize_phased_variants}" == "false" else 'output_dir/phab_bench'
        with open(f'{output_dir}/summary.json') as file:
            json_file = json.load(file)

        df = pd.DataFrame({k: v for k,v in json_file.items() if k != 'gt_matrix'}, index=[0])
        df['Base_Name'] = "~{default="" base_sample_name}"
        df['Comp_Name'] = "~{default="" comp_sample_name}"
        df['Experiment'] = "~{experiment}"
        df['Interval'] = "WholeGenome"

        df.to_csv("TruvariSummary.tsv", sep='\t', index=False)

        CODE

        # Create tables for making closest data and for interval overlap stats
        # Transform VCFs into bed files for more convenient analysis
        # Use POS0/END0 for coherence with bed coordinates later
        echo -e "CHROM\tPOS\tEND\tSVLEN\tSVTYPE\tFILTER\tGTMATCH" > input_header.txt
        INPUT_QUERY="%CHROM\t%POS0\t%END0\t%SVLEN\t%SVTYPE\t%FILTER\n"

        bcftools query -i 'INFO/SVTYPE!="."' -f"${INPUT_QUERY}" base.vcf.gz | bedtools sort -i - > base-preheader.bed
        cat input_header.txt base-preheader.bed > base.bed

        bcftools query -i 'INFO/SVTYPE!="."' -f"${INPUT_QUERY}" comp.vcf.gz | bedtools sort -i - > comp-preheader.bed
        cat input_header.txt comp-preheader.bed > comp.bed

        echo -e "CHROM\tPOS\tEND\tSVLEN\tSVTYPE\tFILTER\tGTMATCH" > truvari_header.txt
        MAIN_QUERY="%CHROM\t%POS0\t%END0\t%SVLEN\t%SVTYPE\t%FILTER\t%GTMatch\n"

        bcftools query -i 'INFO/SVTYPE!="."' -f"${MAIN_QUERY}" "${RESULTS_STEM}tp-base.vcf.gz" | bedtools sort -i - > tp-base-preheader.bed
        cat truvari_header.txt tp-base-preheader.bed > tp-base.bed

        bcftools query -i 'INFO/SVTYPE!="."' -f"${MAIN_QUERY}" "${RESULTS_STEM}fn.vcf.gz" | bedtools sort -i - > fn-preheader.bed
        cat truvari_header.txt fn-preheader.bed > fn.bed

        bcftools query -i 'INFO/SVTYPE!="."' -f"${MAIN_QUERY}" "${RESULTS_STEM}tp-comp.vcf.gz" | bedtools sort -i - > tp-comp-preheader.bed
        cat truvari_header.txt tp-comp-preheader.bed > tp-comp.bed

        bcftools query -i 'INFO/SVTYPE!="."' -f"${MAIN_QUERY}" "${RESULTS_STEM}fp.vcf.gz" | bedtools sort -i - > fp-preheader.bed
        cat truvari_header.txt fp-preheader.bed > fp.bed

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.1"
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        memory: runtimeAttributes.memory + " GB"
        cpu: runtimeAttributes.cpu
    }

    output {
        File fn = "output_dir/fn.vcf.gz"
        File fn_index = "output_dir/fn.vcf.gz.tbi"
        File fp = "output_dir/fp.vcf.gz"
        File fp_index = "output_dir/fp.vcf.gz.tbi"
        File tp_base = "output_dir/tp-base.vcf.gz"
        File tp_base_index = "output_dir/tp-base.vcf.gz.tbi"
        File tp_comp = "output_dir/tp-comp.vcf.gz"
        File tp_comp_index = "output_dir/tp-comp.vcf.gz.tbi"

        File base_table = "base.bed"
        File comp_table = "comp.bed"
        File tp_base_table = "tp-base.bed"
        File fn_table = "fn.bed"
        File tp_comp_table = "tp-comp.bed"
        File fp_table = "fp.bed"

        File summary = "TruvariSummary.tsv"
    }
}

task CollectTruvariClosestStats {
    input {
        File base_table
        File comp_table
        File fp_table
        File fn_table

        String? base_sample_name
        String? comp_sample_name

        String experiment = ""

        File ref_fai

        Int k_closest = 3    # Number of close by variants to compare to

        # Runtime parameters
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(4 * size(base_table, "GB")) + 100, "cpu": 4, "memory": 16}
    }

    command <<<
        set -xueo pipefail

        # Generate "closest" data
        echo -e "BASENAME\tCOMPNAME\tExperiment\tLCHROM\tLPOS\tLEND\tLLEN\tLTYPE\tLFILTER\tLGTMatch\tRCHROM\tRPOS\tREND\tRLEN\tRTYPE\tRFILTER\tDIST" > header.txt

        bedtools closest -a ~{fp_table} -b ~{base_table} -k ~{k_closest} -D ref > fp-closest.bed
        awk -v OFS='\t' '{ print "~{base_sample_name}-Base", "~{comp_sample_name}-Comp", "~{experiment}", $0 }' fp-closest.bed > fp-closest-samples.bed
        cat header.txt fp-closest-samples.bed > fp-closest-final.bed

        bedtools closest -a ~{fn_table} -b ~{comp_table} -k ~{k_closest} -D ref > fn-closest.bed
        awk -v OFS='\t' '{ print "~{base_sample_name}-Base", "~{comp_sample_name}-Comp", "~{experiment}", $0 }' fn-closest.bed > fn-closest-samples.bed
        cat header.txt fn-closest-samples.bed > fn-closest-final.bed

        # Correct POS0/END back to original VCF POS/END coordinates
        python3 << CODE
        import pandas as pd

        for file in ['fp-closest-final.bed', 'fn-closest-final.bed']:
            df = pd.read_csv(file, sep='\t')
            df['LPOS'] = df['LPOS'] + 1
            df['LEND'] = df['LEND'] + 1
            df['RPOS'] = df['RPOS'] + 1
            df['REND'] = df['REND'] + 1

        CODE

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.0"
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        memory: runtimeAttributes.memory + " GB"
        cpu: runtimeAttributes.cpu
    }

    output {
        File fp_closest = "fp-closest-final.bed"
        File fn_closest = "fn-closest-final.bed"
    }
}

task ComputeTruvariIntervalSummaryStats {
    input {
        File tp_base_intervals
        File tp_comp_intervals
        File fp_intervals
        File fn_intervals

        String base_sample_name = ""
        String comp_sample_name = ""

        Array[String] bed_labels

        String experiment = ""

        Array[Int] svlen_bin_cutoffs = [100, 250, 1000, 2500, 10000, 25000]
        Array[Int] overlap_percents = [0, 25, 50, 75, 100]

        # Runtime Parameters
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(4 * size(tp_base_intervals, "GB")) + 100, "cpu": 4, "memory": 16}
    }

    command <<<
        set -xueo pipefail

        python3 << CODE
        import pandas as pd

        tp_base_df = pd.read_csv("~{tp_base_intervals}", sep='\t')
        tp_comp_df = pd.read_csv("~{tp_comp_intervals}", sep='\t')
        fp_df = pd.read_csv("~{fp_intervals}", sep='\t')
        fn_df = pd.read_csv("~{fn_intervals}", sep='\t')

        # Bin counts by SVLEN
        def add_suffix(n):
            if n < 1_000:
                return f'{n}'
            elif n < 1_000_000:
                return f'{n/1000:.1f}k'
            elif n < 1_000_000_000:
                return f'{n/1_000_000:.1f}M'
            else:
                return f'{n}'

        SVLEN_Bin_cutoffs = [~{default="" sep=", " svlen_bin_cutoffs}]
        SVLEN_Bin_start = [f'<{add_suffix(SVLEN_Bin_cutoffs[0])}']
        SVLEN_Bins_middle = [f'{add_suffix(SVLEN_Bin_cutoffs[i])}-{add_suffix(SVLEN_Bin_cutoffs[i+1])}' for i in range(len(SVLEN_Bin_cutoffs)-1)]
        SVLEN_Bin_end = [f'>{add_suffix(SVLEN_Bin_cutoffs[-1])}']
        SVLEN_Bins = SVLEN_Bin_start + SVLEN_Bins_middle + SVLEN_Bin_end

        def bin_svlen(num):
            for i, bin in enumerate(SVLEN_Bin_cutoffs):
                if num < bin:
                    return SVLEN_Bins[i]
            return SVLEN_Bins[-1]

        tp_base_df['SVLEN_Bin'] = tp_base_df['SVLEN'].apply(bin_svlen)
        tp_comp_df['SVLEN_Bin'] = tp_comp_df['SVLEN'].apply(bin_svlen)
        fp_df['SVLEN_Bin'] = fp_df['SVLEN'].apply(bin_svlen)
        fn_df['SVLEN_Bin'] = fn_df['SVLEN'].apply(bin_svlen)

        intervals = ["~{default="" sep="\", \"" bed_labels}"]
        overlap_percents = [~{default="" sep=", " overlap_percents}] + ['>0']
        full_counts_df = pd.DataFrame()
        for interval in intervals:
            for pct in overlap_percents:
                merged_df = pd.DataFrame({'SVTYPE': [], 'SVLEN_Bin': [], 'FILTER': []})
                full_merge_columns = ['SVTYPE', 'SVLEN_Bin', 'FILTER', 'Interval', 'Pct_Overlap']
                for label, df in zip(['tp_base', 'tp_comp', 'fp', 'fn'], [tp_base_df, tp_comp_df, fp_df, fn_df]):
                    for overlap_mode in ['full', 'breakpoint']:
                        for gt_mode in ['match', 'no_match']:
                            if gt_mode == 'no_match' or 'tp' in label:    # Skip match GT for FP/FN
                                if overlap_mode == 'full':
                                    count_label = label
                                    if pct == '>0':
                                        subset_condition = df[f'{interval}-overlap'] > 0
                                    else:
                                        subset_condition = df[f'{interval}-overlap'] >= pct/100
                                else:
                                    count_label = f'{label}-BEND'
                                    if pct == '>0':
                                        subset_condition = (df[f'{interval}-LBEND-overlap'] > 0) | (df[f'{interval}-RBEND-overlap'] > 0)
                                    else:
                                        subset_condition = (df[f'{interval}-LBEND-overlap'] >= pct/100) | (df[f'{interval}-RBEND-overlap'] >= pct/100)

                                if gt_mode == 'match':
                                    subset_condition = subset_condition & (df['GTMATCH'] == 0)

                                sub_df = df[subset_condition].copy()
                                if len(sub_df) > 0:
                                    gt_suffix = '-gt' if gt_mode == 'match' else ''
                                    sub_counts = sub_df.groupby(['SVTYPE', 'SVLEN_Bin', 'FILTER']).apply(len).reset_index().rename(columns={0: f'{count_label}{gt_suffix}_count'})

                                    svtype_all_counts = sub_counts.groupby(['SVLEN_Bin', 'FILTER']).sum(numeric_only=True).reset_index()
                                    svtype_all_counts['SVTYPE'] = 'ALL'
                                    sub_counts = pd.concat([sub_counts, svtype_all_counts])

                                    svlen_all_counts = sub_counts.groupby(['SVTYPE', 'FILTER']).sum(numeric_only=True).reset_index()
                                    svlen_all_counts['SVLEN_Bin'] = 'ALL'
                                    sub_counts = pd.concat([sub_counts, svlen_all_counts])

                                    sub_counts['Pct_Overlap'] = pct
                                    sub_counts['Interval'] = interval

                                    merge_cols = [c for c in full_merge_columns if c in merged_df.columns]
                                    merged_df = merged_df.merge(sub_counts, on=merge_cols, how='outer').fillna(0)
                full_counts_df = pd.concat([full_counts_df, merged_df])

        # Generate Precision, Recall, F1 stats
        full_counts_df['Precision'] = full_counts_df['tp_comp_count'] / (full_counts_df['tp_comp_count'] + full_counts_df['fp_count'])
        full_counts_df['Precision-BEND'] = full_counts_df['tp_comp-BEND_count'] / (full_counts_df['tp_comp-BEND_count'] + full_counts_df['fp-BEND_count'])
        full_counts_df['Recall'] = full_counts_df['tp_base_count'] / (full_counts_df['tp_base_count'] + full_counts_df['fn_count'])
        full_counts_df['Recall-BEND'] = full_counts_df['tp_base-BEND_count'] / (full_counts_df['tp_base-BEND_count'] + full_counts_df['fn-BEND_count'])
        full_counts_df['F1_Score'] = 2 * full_counts_df['Precision'] * full_counts_df['Recall'] / (full_counts_df['Precision'] + full_counts_df['Recall'])
        full_counts_df['F1_Score-BEND'] = 2 * full_counts_df['Precision-BEND'] * full_counts_df['Recall-BEND'] / (full_counts_df['Precision-BEND'] + full_counts_df['Recall-BEND'])

        full_counts_df['GT_concordance'] = full_counts_df['tp_comp-gt_count'] / (full_counts_df['tp_comp_count'])

        # Add sample names
        full_counts_df['Base_Sample_Name'] = "~{base_sample_name}"
        full_counts_df['Comp_Sample_Name'] = "~{comp_sample_name}"

        full_counts_df['Experiment'] = "~{experiment}"

        column_order = ['Base_Sample_Name', 'Comp_Sample_Name', 'Experiment', 'SVTYPE', 'SVLEN_Bin', 'FILTER', 'Interval', 'Pct_Overlap',
                        'tp_base_count', 'tp_base-BEND_count', 'tp_base-gt_count', 'tp_comp_count', 'tp_comp-BEND_count', 'tp_comp-gt_count', 'fp_count', 'fp-BEND_count',
                        'fn_count', 'fn-BEND_count', 'Precision', 'Precision-BEND', 'Recall', 'Recall-BEND', 'F1_Score', 'F1_Score-BEND', 'GT_concordance']

        full_counts_df[column_order].to_csv('FullTruvariSummary.tsv', sep='\t', index=False)

        CODE

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.0"
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        memory: runtimeAttributes.memory + " GB"
        cpu: runtimeAttributes.cpu
    }

    output {
        File full_summary = "FullTruvariSummary.tsv"
    }
}
