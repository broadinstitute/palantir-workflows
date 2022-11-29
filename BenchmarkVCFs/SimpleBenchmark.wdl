version 1.0


# Object holding all reference files, to ensure all localized in tasks
struct Reference {
    File fasta
    File index
    File dict
}

# Object representing expression for subsetting VCFs
struct VariantSelector {
    String bcf_label
    String bcf_selector
}

# Object representing file with intervals to subset analysis over
struct Stratifier {
    String label
    File? interval
}

# Main workflow: performs evaluation of call_vcf against base_vcf using vcfeval comparison engine,
# over different analysis modes.
workflow SimpleBenchmark {
    input {
        # VCF information
        File base_vcf
        File base_vcf_index
        String base_output_sample_name
        String? base_vcf_sample_name

        File call_vcf
        File call_vcf_index
        String call_output_sample_name
        String? call_vcf_sample_name

        # Reference information
        File ref_fasta
        File ref_index
        File ref_dict

        # Subsetting inputs using intervals
        Array[File] strat_intervals = []
        Array[String] strat_labels = []
        Int interval_padding = 0
        String subset_gatk_tag = "4.2.3.0"

        # Subsetting inputs using variant properties; Note SNP & INDEL are separated automatically later
        Array[String] bcf_selectors = []
        Array[String] bcf_labels = []

        # Evaluation inputs
        File? evaluation_bed
        Array[String] score_fields = ["GQ"]

        # Columns to add to output files
        String? experiment_value
        String? extra_column_name
        String? extra_column_value

        # Toggle for more tries on preemptible machines for potentially cheaper runs at the risk of longer runtime
        Int preemptible = 3

        # Null File type -- do NOT assign; Needed until WDL has better support for null File values
        File? NULL_FILE
    }

    # Create reference object
    Reference reference = {"fasta": ref_fasta, "index": ref_index, "dict": ref_dict}

    # Add in default "Whole Genome" empty interval list
    Array[File] full_strat_intervals = flatten([[""], strat_intervals])
    Array[String] full_strat_labels = flatten([[""], strat_labels])

    # Make Stratifier objects
    scatter (strat in zip(full_strat_labels, full_strat_intervals)) {
        Stratifier stratifier_list = {"label": strat.left, "interval": if strat.right != "" then strat.right else NULL_FILE}
    }

    # Make VariantSelector objects, with default None, Het, and HomVar types
    Array[String] default_bcf_selectors =  ["", "GT='het'", "GT='hom'"]
    Array[String] default_bcf_labels = ["", "Het", "HomVar"]
    Array[String] full_bcf_selectors = flatten([default_bcf_selectors, bcf_selectors])
    Array[String] full_bcf_labels = flatten([default_bcf_labels, bcf_labels])
    scatter (selection in zip(full_bcf_labels, full_bcf_selectors)) {
        VariantSelector selector_list = {"bcf_label": selection.left, "bcf_selector": selection.right}
    }

    # Perform analysis over all Stratifiers
    scatter (stratifier in stratifier_list) {
        call SubsetVCF as SubsetEval {
            input:
                input_vcf=call_vcf,
                input_vcf_index=call_vcf_index,
                input_sample_name=call_vcf_sample_name,
                reference=reference,
                stratifier=stratifier,
                interval_padding=interval_padding,
                gatk_tag=subset_gatk_tag,
                preemptible=preemptible
        }

        call SubsetVCF as SubsetTruth {
            input:
                input_vcf=base_vcf,
                input_vcf_index=base_vcf_index,
                input_sample_name=base_vcf_sample_name,
                reference=reference,
                stratifier=stratifier,
                interval_padding=interval_padding,
                gatk_tag=subset_gatk_tag,
                preemptible=preemptible
        }

        # Run over different score_fields to produce other ROC outputs
        scatter (score_field in score_fields) {
            call VCFEval as StandardVCFEval {
                input:
                    call_vcf=SubsetEval.output_vcf,
                    call_vcf_index=SubsetEval.output_vcf_index,
                    call_output_sample_name=call_output_sample_name,
                    base_vcf=SubsetTruth.output_vcf,
                    base_vcf_index=SubsetTruth.output_vcf_index,
                    base_output_sample_name=base_output_sample_name,
                    reference=reference,
                    evaluation_bed=evaluation_bed,
                    score_field=score_field,
                    strat_label=stratifier.label,
                    preemptible=preemptible
            }
        }

        scatter (selector in selector_list) {
            # Use the first outputs for tp_base_vcf, fp_vcf, etc. since these should be the same across all score_fields
            call BCFToolsStats {
                input:
                    tp_base_vcf=StandardVCFEval.tp_base_vcf[0],
                    tp_base_index=StandardVCFEval.tp_base_index[0],
                    tp_call_vcf=StandardVCFEval.tp_call_vcf[0],
                    tp_call_index=StandardVCFEval.tp_call_index[0],
                    fp_vcf=StandardVCFEval.fp_vcf[0],
                    fp_index=StandardVCFEval.fp_index[0],
                    fn_vcf=StandardVCFEval.fn_vcf[0],
                    fn_index=StandardVCFEval.fn_index[0],
                    strat_label=select_first([stratifier.label, ""]),
                    bcf_selector=selector.bcf_selector,
                    bcf_label=selector.bcf_label,
                    call_output_sample_name=call_output_sample_name,
                    base_output_sample_name=base_output_sample_name,
                    evaluation_bed=evaluation_bed
            }
        }
    }

    call CombineSummaries {
        input:
            ROC_summaries=flatten(StandardVCFEval.ROC_summary),
            SN_summaries=flatten(BCFToolsStats.full_sn),
            IDD_summaries=flatten(BCFToolsStats.full_idd),
            ST_summaries=flatten(BCFToolsStats.full_st),
            experiment_value=experiment_value,
            extra_column_name=extra_column_name,
            extra_column_value=extra_column_value

    }

    output {
        File simple_summary = CombineSummaries.simple_summary
        File combined_IDD = CombineSummaries.IDD_combined_summaries
        File combined_ST = CombineSummaries.ST_combined_summaries
        File combined_ROC = CombineSummaries.ROC_combined_summaries
    }

}

task SubsetVCF {
    input {
        File input_vcf
        File input_vcf_index
        String? input_sample_name
        Reference reference

        Stratifier stratifier
        Int interval_padding = 0    # Amount of bases to add around stratification intervals

        String gatk_tag
        Int? preemptible
        Int disk_size = 10 + ceil(4.2 * size(input_vcf, "GB") + 2.2 * size(input_vcf_index, "GB") + size(reference.fasta, "GB"))
        Int cpu = 4
        Int memory = 16
    }

    command <<<
        set -xeuo pipefail

        # Create symlink for VCF & index in case their paths are different, e.g. when using TDR
        ln -s ~{input_vcf} input.vcf.gz
        ln -s ~{input_vcf_index} input.vcf.gz.tbi

        # Subset to given sample / loci; Output subsetted VCF
        # Add variable padding to ensure capture of variants on boundaries downstream
        # Remove unused alternates to make bcftools stats more accurate on multisample callsets
        gatk SelectVariants \
            -V input.vcf.gz \
            -R ~{reference.fasta} \
            ~{"-sn " + input_sample_name} \
            ~{"-L " + stratifier.interval} \
            ~{"--interval-padding " + interval_padding} \
            --remove-unused-alternates \
            -O subset.vcf.gz

    >>>

    runtime{
        docker: "us.gcr.io/broad-gatk/gatk:" + gatk_tag
        preemptible: select_first([preemptible, 0])
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: "16"
        cpu: cpu
        memory: memory + " GB"
    }

    output {
        File output_vcf = "subset.vcf.gz"
        File output_vcf_index = "subset.vcf.gz.tbi"
    }
}

task VCFEval {
    input {
        # Input VCF Files
        File call_vcf
        File call_vcf_index
        String call_output_sample_name
        File base_vcf
        File base_vcf_index
        String base_output_sample_name

        Reference reference

        # Interval File to Subset Analysis on given actual truth data
        File? evaluation_bed

        # String for VCF field to use as ROC score
        String score_field

        # Label to add to ROC table outputs
        String strat_label

        # vcfeval Arguments
        Boolean passing_only = true
        Boolean require_matching_genotypes = true
        Boolean enable_ref_overlap = true

        # Runtime params
        Int? preemptible
        Int disk_size = ceil(size(call_vcf, "GB") + size(base_vcf, "GB") + size(reference.fasta, "GB")) + 25
        Int cpu = 8
        Int memory = 16
        String rtg_docker_version = "v1.0"
    }

    command <<<
        set -xeuo pipefail

        rtg format -o rtg_ref ~{reference.fasta}
        rtg vcfeval \
            ~{false="--all-records" true="" passing_only} \
            ~{false="--squash-ploidy" true="" require_matching_genotypes} \
            ~{true="--ref-overlap" false="" enable_ref_overlap} \
            -b ~{base_vcf} \
            -c ~{call_vcf} \
            ~{"-e " + evaluation_bed} \
            --vcf-score-field="~{score_field}" \
            --output-mode split \
            --decompose \
            --roc-subset snp,mnp,indel \
            -t rtg_ref \
            -o output_dir \

        python3 << CODE
        import gzip
        import pandas as pd

        roc_summary = pd.DataFrame()
        for Type in ['snp', 'mnp', 'indel']:
            file_path = f'output_dir/{Type}_roc.tsv.gz'
            with gzip.open(file_path, 'rt') as file:
                roc = [line for line in file.readlines()]
            header_names = [line.replace('#', '').replace('\n', '') for line in roc if '#' in line][-1].split('\t')
            df = pd.read_csv(file_path, sep='\t', comment='#', header=None, names=header_names)

            rename_columns = {'score': 'Score', 'true_positives_baseline': 'TP_Base',
                  'false_positives': 'FP', 'true_positives_call': 'TP_Call', 'false_negatives': 'FN',
                  'precision': 'Precision', 'sensitivity': 'Recall', 'f_measure': 'F1_Score'}

            df = df.rename(columns=rename_columns)
            df['Type'] = Type.upper()
            df['Stratifier'] = "~{strat_label}"
            df['Score_Field'] = "~{score_field}"
            df['Call_Name'] = "~{call_output_sample_name}"
            df['Base_Name'] = "~{base_output_sample_name}"

            roc_summary = pd.concat([roc_summary, df])

        roc_summary.to_csv('ROC_summary.tsv', sep='\t', index=False)

        CODE

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/rtg:" + rtg_docker_version
        preemptible: select_first([preemptible, 0])
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory + " GB"
    }

    output {
        File ROC_summary = "ROC_summary.tsv"

        File tp_base_vcf = "output_dir/tp-baseline.vcf.gz"
        File tp_base_index = "output_dir/tp-baseline.vcf.gz.tbi"
        File tp_call_vcf = "output_dir/tp.vcf.gz"
        File tp_call_index = "output_dir/tp.vcf.gz.tbi"
        File fp_vcf = "output_dir/fp.vcf.gz"
        File fp_index = "output_dir/fp.vcf.gz.tbi"
        File fn_vcf = "output_dir/fn.vcf.gz"
        File fn_index = "output_dir/fn.vcf.gz.tbi"
    }
}

task BCFToolsStats {
    input {
        # Post-vcfeval files
        File tp_base_vcf
        File tp_base_index
        File tp_call_vcf
        File tp_call_index
        File fp_vcf
        File fp_index
        File fn_vcf
        File fn_index

        String strat_label
        String bcf_selector
        String bcf_label

        String call_output_sample_name
        String base_output_sample_name

        File? evaluation_bed

        Int disk_size = ceil(size(tp_base_vcf, "GB") + size(tp_call_vcf, "GB") + size(fp_vcf, "GB") + size(fn_vcf, "GB")) + 20
        Int cpu = 16
        Int memory = 32
    }

    # Handle parsing bcf_selector with surrounding quotes, with empty case handled separately
    String selection = if bcf_selector!= "" then "-i " + '"' + bcf_selector + '"' else ""

    command <<<
        set -xeuo pipefail

        # Subset to evaluation bed if provided, to ensure not too many FNs picked up outside of it
        bcftools stats ~{selection} -s- ~{"-R" + evaluation_bed} --threads ~{cpu} ~{tp_base_vcf} > tp_base_stats.tsv
        bcftools stats ~{selection} -s- ~{"-R" + evaluation_bed} --threads ~{cpu} ~{tp_call_vcf} > tp_call_stats.tsv
        bcftools stats ~{selection} -s- ~{"-R" + evaluation_bed} --threads ~{cpu} ~{fp_vcf} > fp_stats.tsv
        bcftools stats ~{selection} -s- ~{"-R" + evaluation_bed} --threads ~{cpu} ~{fn_vcf} > fn_stats.tsv

        python3 << CODE
        import io
        import pandas as pd

        def make_dfs(file_path):
            with open(file_path, 'r') as file:
                bcf_file = file.readlines()

            # sample-specific stats
            sn_df = pd.read_csv(io.StringIO('\n'.join([x for x in bcf_file if 'SN\t' in x])), sep='\t')
            final_sn_df = sn_df.drop(columns=['# SN', '[2]id']).rename(columns={'[3]key': 'Category', '[4]value': 'Count'})[3:]
            final_sn_df['Category'] = final_sn_df['Category'].apply(lambda x: x.replace('number of ', '').replace(':', ''))
            idd_df = pd.read_csv(io.StringIO('\n'.join([x for x in bcf_file if 'IDD\t' in x])), sep='\t')
            final_idd_df = idd_df.drop(columns=['# IDD', '[2]id', '[5]number of genotypes', '[6]mean VAF']).rename(columns={
                '[3]length (deletions negative)': 'INDEL_Length', '[4]number of sites': 'Count'
            })
            st_df = pd.read_csv(io.StringIO('\n'.join([x for x in bcf_file if 'ST\t' in x])), sep='\t')
            final_st_df = st_df.drop(columns=['# ST', '[2]id']).rename(columns={'[3]type': 'Substitution', '[4]count': 'Count'})

            # run-specific stats
            qual_df = pd.read_csv(io.StringIO('\n'.join([x for x in bcf_file if 'QUAL\t' in x])), sep='\t')
            final_qual_df = qual_df.drop(columns=['# QUAL', '[2]id']).rename(columns={
                '[3]Quality': 'QUAL', '[4]number of SNPs': 'SNP_Count', '[5]number of transitions (1st ALT)': 'Ti_Count',
                '[6]number of transversions (1st ALT)': 'Tv_Count', '[7]number of indels': 'INDEL_Count'
            })
            dp_df = pd.read_csv(io.StringIO('\n'.join([x for x in bcf_file if 'DP\t' in x])), sep='\t')
            final_dp_df = dp_df.drop(columns=['# DP', '[2]id', '[4]number of genotypes', '[5]fraction of genotypes (%)', '[7]fraction of sites (%)']).rename(columns={
                '[3]bin': 'DP', '[6]number of sites': 'Count'
            })

            return {
                'SN_df': final_sn_df,
                'IDD_df': final_idd_df,
                'ST_df': final_st_df,
                'QUAL_df': final_qual_df,
                'DP_df': final_dp_df
            }

        df_dict = {
            "TP_Base" : make_dfs("tp_base_stats.tsv"),
            "TP_Call": make_dfs("tp_call_stats.tsv"),
            "FP": make_dfs("fp_stats.tsv"),
            "FN": make_dfs("fn_stats.tsv")
        }

        # Make full SN (Summary Numbers) df
        full_sn_df = pd.DataFrame({'Category': df_dict['TP_Base']['SN_df']['Category']})
        for stat in df_dict:
            full_sn_df = df_dict[stat]['SN_df'].rename(columns={'Count': f'{stat}_Count'}).merge(full_sn_df, on='Category', how='outer')
        full_sn_df = full_sn_df.fillna(0)
        full_sn_df['Call_Name'] = "~{call_output_sample_name}"
        full_sn_df['Base_Name'] = "~{base_output_sample_name}"
        full_sn_df['Stratifier'] = "~{strat_label}"
        full_sn_df['BCF_Label'] = "~{bcf_label}"

        # Make full IDD (InDel Distribution) df
        full_idd_df = pd.DataFrame({'INDEL_Length': df_dict['TP_Base']['IDD_df']['INDEL_Length']})
        for stat in df_dict:
            full_idd_df = df_dict[stat]['IDD_df'].rename(columns={'Count': f'{stat}_Count'}).merge(full_idd_df, on='INDEL_Length', how='outer')
        full_idd_df = full_idd_df.fillna(0)
        full_idd_df['Call_Name'] = "~{call_output_sample_name}"
        full_idd_df['Base_Name'] = "~{base_output_sample_name}"
        full_idd_df['Stratifier'] = "~{strat_label}"
        full_idd_df['BCF_Label'] = "~{bcf_label}"

        # Make full ST (Substitution) df
        full_st_df = pd.DataFrame({'Substitution': df_dict['TP_Base']['ST_df']['Substitution']})
        for stat in df_dict:
            full_st_df = df_dict[stat]['ST_df'].rename(columns={'Count': f'{stat}_Count'}).merge(full_st_df, on='Substitution', how='outer')
        full_st_df = full_st_df.fillna(0)
        full_st_df['Call_Name'] = "~{call_output_sample_name}"
        full_st_df['Base_Name'] = "~{base_output_sample_name}"
        full_st_df['Stratifier'] = "~{strat_label}"
        full_st_df['BCF_Label'] = "~{bcf_label}"


        # Write files
        full_sn_df.to_csv('Full_SN.tsv', sep='\t', index=False)
        full_idd_df.to_csv('Full_IDD.tsv', sep='\t', index=False)
        full_st_df.to_csv('Full_ST.tsv', sep='\t', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.0"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory + "GB"
    }

    output {
        File full_sn = "Full_SN.tsv"
        File full_idd = "Full_IDD.tsv"
        File full_st = "Full_ST.tsv"
    }
}

task CombineSummaries {
    input {
        Array[File] ROC_summaries

        Array[File] SN_summaries
        Array[File] IDD_summaries
        Array[File] ST_summaries

        String? experiment_value
        String? extra_column_name
        String? extra_column_value

        Int disk_size = 50
        Int cpu = 2
        Int memory = 8
    }

    command <<<
        set -xeuo pipefail

        python << CODE
        import pandas as pd

        # Concat all ROC summaries into one file
        full_ROC = pd.DataFrame()
        for file in ["~{default="" sep="\", \"" ROC_summaries}"]:
            df = pd.read_csv(file, sep='\t')
            full_ROC = pd.concat([full_ROC, df])

        if "~{experiment_value}" != "":
            full_ROC['Experiment'] = "~{experiment_value}"

        if ("~{extra_column_name}" != "") and ("~{extra_column_value}" != ""):
            full_ROC["~{extra_column_name}"] = "~{extra_column_value}"

        full_ROC.to_csv('Full_ROC.tsv', sep='\t', index=False)

        # Gather all tables for bcftools stats outputs
        full_SN = pd.DataFrame()
        for file in ["~{default="" sep="\", \"" SN_summaries}"]:
            df = pd.read_csv(file, sep='\t')
            full_SN = pd.concat([full_SN, df])

        full_IDD = pd.DataFrame()
        for file in ["~{default="" sep="\", \"" IDD_summaries}"]:
            df = pd.read_csv(file, sep='\t')
            full_IDD = pd.concat([full_IDD, df])

        full_ST = pd.DataFrame()
        for file in ["~{default="" sep="\", \"" ST_summaries}"]:
            df = pd.read_csv(file, sep='\t')
            full_ST = pd.concat([full_ST, df])

        # Label INDEL types
        full_IDD['INDEL_Type'] = full_IDD['INDEL_Length'].apply(lambda x: 'Ins' if x > 0 else 'Del')
        full_IDD['Type'] = full_IDD['BCF_Label'].fillna("") + 'INDEL'
        full_IDD = full_IDD.drop(columns=['BCF_Label'])

        # Label SNP sub types
        def ti_tv(sub):
            if (set(sub.split('>')) == {'A', 'G'}) or (set(sub.split('>')) == {'C', 'T'}):
                return 'Ti'
            else:
                return 'Tv'

        full_ST['Substitution_Type'] = full_ST['Substitution'].apply(ti_tv)
        full_ST['Type'] = full_ST['BCF_Label'].fillna("") + 'SNP'
        full_ST = full_ST.drop(columns=['BCF_Label'])

        # Clean up for simple summary
        simple_renaming = {'SNPs': 'SNP', 'MNPs': 'MNP', 'indels': 'INDEL', 'others': 'Other', 'multiallelic sites': 'MA', 'multiallelic SNP sites': 'MASNP'}
        simple_summary = full_SN.copy()
        simple_summary['Category'] = simple_summary['Category'].replace(simple_renaming)
        simple_summary['Type'] = simple_summary['BCF_Label'].fillna("") + simple_summary['Category']
        simple_summary = simple_summary.drop(columns=['BCF_Label', 'Category'])

        # Compute simple stats
        for df in [full_SN, full_IDD, full_ST, simple_summary]:
            df.columns = [c.replace('_Count', '') for c in df.columns]
            df['Precision'] = df['TP_Call'] / (df['TP_Call'] + df['FP'])
            df['Recall'] = df['TP_Base'] / (df['TP_Base'] + df['FN'])
            df['F1_Score'] = 2 * df['Precision'] * df['Recall'] / (df['Precision'] + df['Recall'])

        # Add optional labels
        for df in [full_SN, full_IDD, full_ST, simple_summary]:
            if "~{experiment_value}" != "":
                df['Experiment'] = "~{experiment_value}"

            if ("~{extra_column_name}" != "") and ("~{extra_column_value}" != ""):
                df["~{extra_column_name}"] = "~{extra_column_value}"

        # Reorder columns
        metadata_cols = ['Call_Name', 'Base_Name', 'Stratifier', 'Type']
        metadata_cols = metadata_cols + ["~{extra_column_name}"] if "~{extra_column_name}" != "" else metadata_cols
        metadata_cols = ['Experiment'] + metadata_cols if "~{experiment_value}" != "" else metadata_cols
        stat_cols = ['TP_Call', 'TP_Base', 'FP', 'FN', 'Precision', 'Recall', 'F1_Score']

        full_IDD = full_IDD[metadata_cols + ['INDEL_Type', 'INDEL_Length'] + stat_cols]
        full_ST = full_ST[metadata_cols + ['Substitution', 'Substitution_Type'] + stat_cols]
        simple_summary = simple_summary[metadata_cols + stat_cols]

        # Skip returning full_SN since it is redundant with simple_summary
        simple_summary.to_csv('SimpleSummary.tsv', sep='\t', index=False)
        full_IDD.to_csv('Full_IDD.tsv', sep='\t', index=False)
        full_ST.to_csv('Full_ST.tsv', sep='\t', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory + "GB"
    }

    output {
        File simple_summary = "SimpleSummary.tsv"
        File IDD_combined_summaries = "Full_IDD.tsv"
        File ST_combined_summaries = "Full_ST.tsv"
        File ROC_combined_summaries = "Full_ROC.tsv"
    }
}
