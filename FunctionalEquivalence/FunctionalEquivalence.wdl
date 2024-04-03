version 1.0

#import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/BenchmarkVCFs/SimpleBenchmark.wdl" as BenchmarkVCFs
import "../BenchmarkVCFs/SimpleBenchmark.wdl" as BenchmarkVCFs
import "subworkflows/F1Evaluation.wdl" as F1Evaluation
import "subworkflows/PlotROC.wdl" as PlotROC

struct VcfFile {
    String sample_id
    File file
    File index
    String dataset
    File confidence_intervals
    Int num    # placeholder for position in list for cross product
}

struct TruthVcf {
    String sample_id
    File file
    File index
    File confidence_intervals
    String dataset
}

workflow FunctionalEquivalence {
    input {
        Array[String] sample_id
        Array[String] dataset

        Array[String] confidence_intervals

        Array[String] tool1_vcf
        Array[String] tool1_vcf_index
        Array[String] tool2_vcf
        Array[String] tool2_vcf_index

        Array[String] truth_vcf
        Array[String] truth_vcf_index
        Array[String] truth_vcf_sample_names

        Array[File]? stratifier_intervals
        Array[String]? stratifier_labels

        File ref_fasta
        File ref_index
        File haplotype_map

        String tool1_label
        String tool2_label
        String? additional_label

        Boolean signed_difference = false

        Boolean passingOnly = true
        Boolean requireMatchingGenotypes = true
        String vcfScoreField = "QUAL"
        Int? threadsVcfEval = 2
        Int? preemptible = 3
    }

    scatter (i in range(length(tool1_vcf))) {
        VcfFile tool1_inputs = {"sample_id": sample_id[i], "file": tool1_vcf[i], "index": tool1_vcf_index[i], "dataset": dataset[i], "confidence_intervals": confidence_intervals[i], "num": i}
    }

    scatter (i in range(length(tool2_vcf))) {
        VcfFile tool2_inputs = {"sample_id": sample_id[i], "file": tool2_vcf[i], "index": tool2_vcf_index[i], "dataset": dataset[i], "confidence_intervals": confidence_intervals[i], "num": i}
    }

    scatter (i in range(length(truth_vcf))) {
        TruthVcf truth_inputs = {"sample_id": truth_vcf_sample_names[i], "file": truth_vcf[i], "index": truth_vcf_index[i], "dataset": dataset[i], "confidence_intervals": confidence_intervals[i]}
    }

    ## Evaluate against the truth files
    # Only used for F1Evaluation half of pipeline (ROC data)
    scatter (paired_vcfs in zip(tool1_inputs, truth_inputs)) {
        call BenchmarkVCFs.SimpleBenchmark as EvalVsTruthTool1 {
            input:
                base_vcf=paired_vcfs.right.file,
                base_vcf_index=paired_vcfs.right.index,
                base_output_sample_name=paired_vcfs.right.sample_id,
                base_vcf_sample_name=paired_vcfs.right.sample_id,
                query_vcf=paired_vcfs.left.file,
                query_vcf_index=paired_vcfs.left.index,
                query_output_sample_name=paired_vcfs.left.sample_id,
                query_vcf_sample_name=paired_vcfs.left.sample_id,
                ref_fasta=ref_fasta,
                ref_index=ref_index,
                haplotype_map=haplotype_map,
                stratifier_intervals=stratifier_intervals,
                stratifier_labels=stratifier_labels,
                evaluation_intervals=paired_vcfs.right.confidence_intervals,
                score_field="QUAL",
                experiment="EvalVsTruthTool1",
                extra_column_names=["Dataset", "Replicate", "Tool", "Interval-test"],
                extra_column_values=[paired_vcfs.left.dataset, paired_vcfs.left.num, tool1_label, "WholeGenome"],
        }
    }

    scatter (paired_vcfs in zip(tool2_inputs, truth_inputs)) {
        call BenchmarkVCFs.SimpleBenchmark as EvalVsTruthTool2 {
            input:
                base_vcf=paired_vcfs.right.file,
                base_vcf_index=paired_vcfs.right.index,
                base_output_sample_name=paired_vcfs.right.sample_id,
                base_vcf_sample_name=paired_vcfs.right.sample_id,
                query_vcf=paired_vcfs.left.file,
                query_vcf_index=paired_vcfs.left.index,
                query_output_sample_name=paired_vcfs.left.sample_id,
                query_vcf_sample_name=paired_vcfs.left.sample_id,
                ref_fasta=ref_fasta,
                ref_index=ref_index,
                haplotype_map=haplotype_map,
                stratifier_intervals=stratifier_intervals,
                stratifier_labels=stratifier_labels,
                evaluation_intervals=paired_vcfs.right.confidence_intervals,
                score_field="QUAL",
                experiment="EvalVsTruthTool2",
                extra_column_names=["Dataset", "Replicate", "Tool", "Interval-test"],
                extra_column_values=[paired_vcfs.left.dataset, paired_vcfs.left.num, tool2_label, "WholeGenome"]
        }
    }

    # Pair sample IDs with pairs of ROC table outputs from above two benchmarks
    Array[Pair[String, Pair[File, File]]] plot_roc_tables = zip(sample_id, zip(EvalVsTruthTool1.ROCStats, EvalVsTruthTool2.ROCStats))

    # Make ROC plot for each sample ID with paired ROC table outputs between the two tools
    scatter(table in plot_roc_tables) {
        call PlotROC.PlotROC as PlotROC {
            input:
                sample_id = table.left,
                roc_tables = [table.right.left, table.right.right],
                tool1_label = tool1_label,
                tool2_label = tool2_label,
                additional_label = additional_label,
                preemptible = preemptible
        }
    }

    ## Evaluate across the two tools
    # Only used for FEEvaluation half of pipeline
    scatter (paired_vcfs in zip(tool1_inputs, tool2_inputs)) {
        call BenchmarkVCFs.SimpleBenchmark as EvalInterTool {
            input:
                base_vcf=paired_vcfs.right.file,
                base_vcf_index=paired_vcfs.right.index,
                base_output_sample_name=paired_vcfs.right.sample_id,
                base_vcf_sample_name=paired_vcfs.right.sample_id,
                query_vcf=paired_vcfs.left.file,
                query_vcf_index=paired_vcfs.left.index,
                query_output_sample_name=paired_vcfs.left.sample_id,
                query_vcf_sample_name=paired_vcfs.left.sample_id,
                ref_fasta=ref_fasta,
                ref_index=ref_index,
                haplotype_map=haplotype_map,
                evaluation_intervals=paired_vcfs.left.confidence_intervals,
                score_field="QUAL",
                stratifier_intervals=stratifier_intervals,
                stratifier_labels=stratifier_labels,
                experiment="EvalInterTool",
                extra_column_names=["Dataset", "Replicate"],
                extra_column_values=[paired_vcfs.left.dataset, paired_vcfs.left.num]
        }
    }

    ## Evaluate within the same tool all possible pairs for both tools
    scatter (index in cross(range(length(tool1_inputs)), range(length(tool1_inputs)))) {
        if (index.left < index.right) {    # Only check when first has index less than second in cross product so no repeats
            call BenchmarkVCFs.SimpleBenchmark as EvalIntraTool1 {
                input:
                    base_vcf=tool1_inputs[index.left].file,
                    base_vcf_index=tool1_inputs[index.left].index,
                    base_output_sample_name=tool1_inputs[index.left].sample_id,
                    base_vcf_sample_name=tool1_inputs[index.left].sample_id,
                    query_vcf=tool1_inputs[index.right].file,
                    query_vcf_index=tool1_inputs[index.right].index,
                    query_output_sample_name=tool1_inputs[index.right].sample_id,
                    query_vcf_sample_name=tool1_inputs[index.right].sample_id,
                    ref_fasta=ref_fasta,
                    ref_index=ref_index,
                    haplotype_map=haplotype_map,
                    evaluation_intervals=tool1_inputs[index.left].confidence_intervals,
                    score_field="QUAL",
                    stratifier_intervals=stratifier_intervals,
                    stratifier_labels=stratifier_labels,
                    experiment="EvalIntraTool1",
                    extra_column_names=["Dataset", "Replicate"],
                    extra_column_values=[tool1_inputs[index.left].dataset, tool1_inputs[index.left].num]
            }
        }
    }

    scatter (index in cross(range(length(tool2_inputs)), range(length(tool2_inputs)))) {
        if (index.left < index.right) {    # Only check when first has index less than second in cross product so no repeats
            call BenchmarkVCFs.SimpleBenchmark as EvalIntraTool2 {
                input:
                    base_vcf=tool2_inputs[index.left].file,
                    base_vcf_index=tool2_inputs[index.left].index,
                    base_output_sample_name=tool2_inputs[index.left].sample_id,
                    base_vcf_sample_name=tool2_inputs[index.left].sample_id,
                    query_vcf=tool2_inputs[index.right].file,
                    query_vcf_index=tool2_inputs[index.right].index,
                    query_output_sample_name=tool2_inputs[index.right].sample_id,
                    query_vcf_sample_name=tool2_inputs[index.right].sample_id,
                    ref_fasta=ref_fasta,
                    ref_index=ref_index,
                    haplotype_map=haplotype_map,
                    evaluation_intervals=tool2_inputs[index.left].confidence_intervals,
                    score_field="QUAL",
                    stratifier_intervals=stratifier_intervals,
                    stratifier_labels=stratifier_labels,
                    experiment="EvalIntraTool2",
                    extra_column_names=["Dataset", "Replicate"],
                    extra_column_values=[tool2_inputs[index.left].dataset, tool2_inputs[index.left].num]
            }
        }
    }

    Array[File] fe_eval_summaries = flatten([select_all(EvalInterTool.SimpleSummary), select_all(EvalIntraTool1.SimpleSummary), select_all(EvalIntraTool2.SimpleSummary)])

    call FEEvaluation {
        input:
            benchmark_summaries=fe_eval_summaries,
            tool1_label=tool1_label,
            tool2_label=tool2_label,
            additional_label=additional_label
    }

    Array[File] f1_roc_tables = flatten([EvalVsTruthTool1.ROCStats, EvalVsTruthTool2.ROCStats])

    call F1Evaluation.F1Evaluation as F1Evaluation {
        input:
            roc_tables=f1_roc_tables,
            tool1_label=tool1_label,
            tool2_label=tool2_label,
            additional_label=additional_label
    }

    Int fe_status_combined = if FEEvaluation.fe_status > F1Evaluation.fe_status then FEEvaluation.fe_status else F1Evaluation.fe_status


    # Also combine all plots into one image
    call MergePNGs as MergeFE {
        input:
            pngs=FEEvaluation.fe_plots,
            preemptible=2
    }

    call MergePNGs as MergeF1 {
        input:
            pngs=F1Evaluation.f1_plots,
            preemptible=preemptible
    }

    call MergePNGs as MergeROC {
        input:
            pngs = flatten(PlotROC.plots),
            preemptible = preemptible
    }

    call CreateHTMLReport {
        input:
            merged_fe_plots=MergeFE.plots,
            merged_f1_plots=MergeF1.plots,
            fe_status=fe_status_combined,
            additional_label=additional_label,
            preemptible=preemptible
    }


    output {
        Array[File] fe_plots = FEEvaluation.fe_plots
        Array[File] f1_plots = F1Evaluation.f1_plots
        Array[File] roc_plots = flatten(PlotROC.plots)
        File merged_fe_plots = MergeFE.plots
        File merged_f1_plots = MergeF1.plots
        File merged_roc_plots = MergeROC.plots
        File fe_summary = FEEvaluation.fe_summary
        File f1_summary = F1Evaluation.f1_summary
        Int fe_status = fe_status_combined
        File html_report = CreateHTMLReport.report
    }
}


task MergePNGs {
    input {
        Array[File] pngs
        Int? preemptible
    }

    command {
        convert -append ~{sep=" " pngs} plots.png
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/functionalequivalence/merge_pngs:1.0.0"
        preemptible: select_first([preemptible, 0])
        memory: "8 GB"
        disks: "local-disk 20 HDD"
    }

    output{
        File plots = "plots.png"
    }
}

task FEEvaluation {
    input {
        Array[File] benchmark_summaries
        String tool1_label
        String tool2_label
        String? additional_label
    }

    String title_label = if (defined(additional_label)) then additional_label else ", "

    command <<<
        set -xueo pipefail

        python3 << CODE
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import matplotlib

        ## A few settings for the rest of the script
        VARIANT_TYPES = ['SNP', 'INDEL']
        matplotlib.rcParams['text.usetex'] = False
        matplotlib.rcParams['mathtext.default'] = 'regular'
        matplotlib.rcParams['font.family'] = 'serif'

        full_df = pd.DataFrame()
        for file in ["~{sep="\", \"" benchmark_summaries}"]:
            df = pd.read_csv(file, sep="\t")
            full_df = pd.concat([full_df, df])

        ## Collect interval names and sort so WholeGenome is always first
        intervals = full_df['Interval'].unique()
        intervals = ['WholeGenome'] + sorted([i for i in intervals if i != 'WholeGenome'])
        full_df['Interval'] = pd.Categorical(full_df['Interval'], intervals)
        full_df = full_df.sort_values('Interval')

        ## Subset/augment the full_df a bit
        full_df = full_df[full_df['Type'].isin(VARIANT_TYPES)]
        full_df['Jaccard'] = full_df['TP_Base'] / (full_df['TP_Base'] + full_df['FP'] + full_df['FN'])

        ## Create the plot_df from stats of Jaccard
        mean_df = full_df.groupby(['Interval', 'Type', 'Experiment', 'Dataset'])['Jaccard'].mean().reset_index().rename(columns={'Jaccard': 'Jaccard_mean'})
        print(f'Means of Jaccard are: {mean_df}')
        sem_df = full_df.groupby(['Interval', 'Type', 'Experiment', 'Dataset'])['Jaccard'].sem().reset_index().rename(columns={'Jaccard': 'Jaccard_sem'})
        print(f'Stderr of Jaccard are: {sem_df}')
        plot_df = mean_df.merge(sem_df, on=['Interval', 'Type', 'Experiment', 'Dataset'])

        ## Function to fill in FE plot for given parameters and return fe_status indicator value
        def make_grid_plot(ax, plot_df, var_type, stratifier):
            fe_status = 0
            x = np.array([0, 1, 2])
            xticklabels = ["~{tool1_label}", 'Inter', "~{tool2_label}"]

            sub_df = plot_df[(plot_df['Type'] == var_type) & (plot_df['Interval'] == stratifier)]
            sub_df['Experiment'] = pd.Categorical(sub_df['Experiment'], ['EvalIntraTool1', 'EvalInterTool', 'EvalIntraTool2'])
            sub_df = sub_df.sort_values(by='Experiment')
            print("sub_df Experiment sort DEBUG!")
            print(sub_df)

            y = sub_df['Jaccard_mean']
            y_err = sub_df['Jaccard_sem']

            ax.set_xticks(x)
            ax.set_xticklabels(xticklabels)
            ax.set_ylabel('Jaccard score')
            ax.errorbar(x, y, yerr=y_err, c='k', fmt='o', capsize=15)
            ax.set_xlim(-0.5, 2.5)
            ax.grid(axis='y')

            inter_mean = sub_df[sub_df['Experiment'] == "EvalInterTool"]['Jaccard_mean'].values[0]
            inter_sem = sub_df[sub_df['Experiment'] == "EvalInterTool"]['Jaccard_sem'].values[0]
            tool1_mean = sub_df[sub_df['Experiment'] == "EvalIntraTool1"]['Jaccard_mean'].values[0]
            tool1_sem = sub_df[sub_df['Experiment'] == "EvalIntraTool1"]['Jaccard_sem'].values[0]
            tool2_mean = sub_df[sub_df['Experiment'] == "EvalIntraTool2"]['Jaccard_mean'].values[0]
            tool2_sem = sub_df[sub_df['Experiment'] == "EvalIntraTool2"]['Jaccard_sem'].values[0]

            if inter_mean <= max(tool1_mean, tool2_mean):
                titlecolor = 'orange'
                fe_status = 2
            elif inter_mean - np.nan_to_num(inter_sem) < max(tool1_mean + np.nan_to_num(tool1_sem), tool2_mean + np.nan_to_num(tool2_sem)):
                titlecolor = 'yellow'
                fe_status = 1
            else:
                titlecolor = 'white'
            ax.set_title(f'{var_type} {stratifier}', backgroundcolor=titlecolor, zorder=0)
            return fe_status


        ## Make plots
        for dataset in plot_df['Dataset'].unique():
            num_columns = max(3, len(intervals))   # Why cap at 3?
            fig, axes = plt.subplots(len(VARIANT_TYPES), num_columns, figsize=(3*num_columns , 6))

            fe_status = 0
            dataset_plot_df = plot_df[plot_df['Dataset'] == dataset]
            for row, var_type in enumerate(VARIANT_TYPES):
                for col, stratifier in enumerate(intervals[:3]):
                    fe_status = max(fe_status, make_grid_plot(axes[row, col], dataset_plot_df, var_type, stratifier))

            tool1_conc_count = ""
            tool2_conc_count = ""
            inter_conc_count = ""
            concordance_text = f'Concordance (# values ~{tool1_label}: {tool1_conc_count} / ~{tool2_label}: {tool2_conc_count} / inter: {inter_conc_count}'
            title = f'Dataset: {dataset} ~{title_label} ~{tool1_label} vs ~{tool2_label}\n{concordance_text}'
            fig.suptitle(title)
            fig.tight_layout()
            fig.savefig(f'fe_plot_{dataset}.png', dpi=100)

        with open('fe_status.txt', 'w') as file:
            file.write(f'{fe_status}')

        full_df.to_csv('fe_summary.tsv', sep='\t', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim-plots:1.0"
        preemptible: select_first([2, 0])
        memory: 8 + " GB"
        disks: "local-disk 20 HDD"
    }

    output {
        Array[File] fe_plots = glob("fe_plot_*.png")
        File fe_summary = "fe_summary.tsv"
        Int fe_status = read_int("fe_status.txt")
    }
}

task CreateHTMLReport {
    input {
        File merged_fe_plots
        File merged_f1_plots
        Int fe_status
        String additional_label = ""
        Int? preemptible
    }

    String fe_status_string = if fe_status == 0 then "The analysis suggests that the tools are functionally equivalent." else "The analysis suggests that the tools are NOT functionally equivalent."
    String fe_status_color = if fe_status == 0 then "#008000" else "#FF0000"

    command <<<
        set -xeuo pipefail

        fe_plots_base64=$(base64 -w 0 ~{merged_fe_plots})
        f1_plots_base64=$(base64 -w 0 ~{merged_f1_plots})

        cat <<EOF > report.html
<!DOCTYPE html>
<html>
    <head>
        <meta charset="UTF-8">
        <title>FE Report ~{additional_label}</title>
        <style>
            body {
                font-family: sans-serif;
                font-size: 14px;
                padding: 0 26px;
                line-height: 1.6;
            }
            img {
                max-width: 100%;
                max-height: 100%;
            }
            table {
                border-collapse: collapse;
            }
            th, td {
                padding: 5px 10px;
            }
            table, td {
                border: 1px solid black;
            }
        </style>
    </head>
    <body>
        <h2 style="color: ~{fe_status_color};">~{fe_status_string}</h2>
        <h3>For more information about how to interpret the plots, please refer to the <a href="https://github.com/broadinstitute/palantir-workflows/tree/main/FunctionalEquivalence">documentation on GitHub</a>.</h3>
        <h2>FE plots</h2>
        <table>
            <tr>
                <th style="text-align: center;">~{additional_label}</th>
            </tr>
            <tr>
                <td><img src="data:image/png;base64,$fe_plots_base64" /></td>
            </tr>
        </table>
        <h2>F1 plots</h2>
        <table>
            <tr>
                <th style="text-align: center;">~{additional_label}</th>
            </tr>
            <tr>
                <td><img src="data:image/png;base64,$f1_plots_base64" /></td>
            </tr>
        </table>
    </body>
</html>
EOF
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim-plots:1.0"
        preemptible: select_first([preemptible, 0])
        memory: "2 GB"
        disks: "local-disk 20 HDD"
    }

    output {
        File report = "report.html"
    }
}