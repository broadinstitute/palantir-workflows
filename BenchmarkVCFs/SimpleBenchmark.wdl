version 1.0

import "../Utilities/WDLs/IntervalList2Bed.wdl" as IntervalList2Bed
import "../Utilities/WDLs/CreateIGVSession.wdl" as IGV
import "../Utilities/WDLs/MatchFingerprints.wdl" as Fingerprint

# Object holding configuration for runtime parameters to shorten number of optional inputs
struct RuntimeAttributes {
    Int disk_size
    Int cpu
    Int memory
}

# Object holding all reference files, to ensure all localized in tasks
struct Reference {
    File fasta
    File index
}

# Object representing expression for subsetting VCFs into genotype categories
struct GenotypeSelector {
    String bcf_genotype_label
    String bcf_genotype
}

# Object representing file with intervals to subset analysis over
struct StratifierInterval {
    String label
    File intervals
}

# Main workflow: performs evaluation of query_vcf against base_vcf using vcfeval comparison engine,
# over different analysis modes.
workflow SimpleBenchmark {
    input {
        # VCF information
        File base_vcf
        File base_vcf_index
        String base_output_sample_name
        String? base_vcf_sample_name

        File query_vcf
        File query_vcf_index
        String query_output_sample_name
        String? query_vcf_sample_name

        # Reference information
        File ref_fasta
        File ref_index

        # Subsetting inputs using intervals
        Array[File] stratifier_intervals = []
        Array[String] stratifier_labels = []

        # Evaluation inputs
        File? evaluation_intervals
        String score_field = "GQ"

        # Columns to add to output files
        String? experiment
        String? extra_column_name
        String? extra_column_value

        Boolean check_fingerprint = true
        File? haplotype_map

        Boolean create_igv_session = false
        Array[File]? optional_igv_bams
        String igv_session_name = "igv_session"

        # Toggle for more tries on preemptible machines for potentially cheaper runs at the risk of longer runtime
        Int preemptible = 3
    }

    # Create reference object
    Reference reference = {"fasta": ref_fasta, "index": ref_index}

    # Convert any provided interval_lists to beds
    if (defined(evaluation_intervals)) {
        call IntervalList2Bed.IntervalList2Bed as ConvertEvalIntervals {
            input:
                interval_files=select_all([evaluation_intervals]),
                interval_labels=["Evaluation"]
        }

        File converted_evaluation_bed = select_first(ConvertEvalIntervals.bed_files)
    }

    call IntervalList2Bed.IntervalList2Bed as ConvertIntervals {
        input:
            interval_files=stratifier_intervals,
            interval_labels=stratifier_labels
    }

    # Make StratifierInterval objects
    scatter (interval_pair in zip(ConvertIntervals.bed_labels, ConvertIntervals.bed_files)) {
        StratifierInterval stratifier_list = {"label": interval_pair.left, "intervals": interval_pair.right}
    }

    # Collect genotype categories to stratify by in bcftools stats
    Array[String] bcf_genotypes = ["", "GT=\"het\"", "GT=\"AA\""]
    Array[String] bcf_genotype_labels = ["", "Het", "HomVar"]

    scatter (selection in zip(bcf_genotype_labels, bcf_genotypes)) {
        GenotypeSelector genotype_selector_list = {"bcf_genotype_label": selection.left, "bcf_genotype": selection.right}
    }

    if (check_fingerprint) {
        call Fingerprint.MatchFingerprints as CheckFingerprint {
            input:
                input_files=[query_vcf],
                input_indices=[query_vcf_index],
                second_input_files=[base_vcf],
                second_input_indices=[base_vcf_index],
                haplotype_map=select_first([haplotype_map]),
                fail_on_mismatch=true
        }
    }

    call VCFEval as StandardVCFEval {
        input:
            query_vcf=query_vcf,
            query_vcf_index=query_vcf_index,
            query_output_sample_name=query_output_sample_name,
            query_vcf_sample_name=query_vcf_sample_name,
            base_vcf=base_vcf,
            base_vcf_index=base_vcf_index,
            base_output_sample_name=base_output_sample_name,
            base_vcf_sample_name=base_vcf_sample_name,
            reference=reference,
            evaluation_bed=converted_evaluation_bed,
            score_field=score_field,
            preemptible=preemptible
    }

    scatter (genotype_selector in genotype_selector_list) {
        call BCFToolsStats as WholeGenomeStats {
            input:
                combined_vcfeval_output=StandardVCFEval.combined_output,
                combined_vcfeval_output_index=StandardVCFEval.combined_output_index,
                stratifier_label="WholeGenome",
                bcf_genotype=genotype_selector.bcf_genotype,
                bcf_genotype_label=genotype_selector.bcf_genotype_label,
                query_output_sample_name=query_output_sample_name,
                base_output_sample_name=base_output_sample_name
        }
    }

    scatter (subset_condition in cross(stratifier_list, genotype_selector_list)) {
        StratifierInterval stratifier = subset_condition.left
        GenotypeSelector genotype_selector = subset_condition.right
        call BCFToolsStats as SubsetStats {
            input:
                combined_vcfeval_output=StandardVCFEval.combined_output,
                combined_vcfeval_output_index=StandardVCFEval.combined_output_index,
                stratifier_interval=stratifier.intervals,
                stratifier_label=stratifier.label,
                bcf_genotype=genotype_selector.bcf_genotype,
                bcf_genotype_label=genotype_selector.bcf_genotype_label,
                query_output_sample_name=query_output_sample_name,
                base_output_sample_name=base_output_sample_name
        }
    }

    call CombineSummaries {
        input:
            ROC_summaries=[StandardVCFEval.ROC_summary],
            SN_summaries=select_all(flatten([WholeGenomeStats.full_sn, SubsetStats.full_sn])),
            IDD_summaries=select_all(flatten([WholeGenomeStats.full_idd, SubsetStats.full_idd])),
            ST_summaries=select_all(flatten([WholeGenomeStats.full_st, SubsetStats.full_st])),
            experiment=experiment,
            extra_column_name=extra_column_name,
            extra_column_value=extra_column_value
    }

    if (create_igv_session) {
        call IGV.CreateIGVSession as IGVSession {
            input:
            bams=optional_igv_bams,
            vcfs=[StandardVCFEval.combined_output],
            interval_files=stratifier_intervals,
            reference=ref_fasta,
            output_name=igv_session_name
        }
    }


    output {
        File SimpleSummary = CombineSummaries.simple_summary
        File IndelDistributionStats = CombineSummaries.IDD_combined_summaries
        File SNPSubstitutionStats = CombineSummaries.ST_combined_summaries
        File ROCStats = CombineSummaries.ROC_combined_summaries

        File? igv_session = IGVSession.igv_session
    }

}

task VCFEval {
    input {
        # Input VCF Files
        File query_vcf
        File query_vcf_index
        String query_output_sample_name
        String? query_vcf_sample_name
        File base_vcf
        File base_vcf_index
        String base_output_sample_name
        String? base_vcf_sample_name

        Reference reference

        # Interval File to Subset Analysis on given actual truth data
        File? evaluation_bed

        # Par File
        File? par_bed

        # String for VCF field to use as ROC score
        String score_field

        # vcfeval Arguments
        Boolean passing_only = true
        Boolean require_matching_genotypes = true
        Boolean enable_ref_overlap = false

        # Runtime params
        Int? preemptible
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(query_vcf, "GB") + 2 * size(base_vcf, "GB") + size(reference.fasta, "GB")) + 50,
                                                  "cpu": 8, "memory": 16}
    }

    command <<<
        set -xeuo pipefail

        # Normal situation without any PAR bed file
        if [ -z ~{par_bed} ];
        then
            rtg format -o rtg_ref ~{reference.fasta}
            rtg vcfeval \
                ~{false="--all-records" true="" passing_only} \
                ~{false="--squash-ploidy" true="" require_matching_genotypes} \
                ~{true="--ref-overlap" false="" enable_ref_overlap} \
                -b ~{base_vcf} \
                -c ~{query_vcf} \
                ~{"-e " + evaluation_bed} \
                --vcf-score-field="~{score_field}" \
                --output-mode combine \
                --decompose \
                --roc-subset snp,indel \
                -t rtg_ref \
                ~{"--sample " + base_vcf_sample_name + "," + query_vcf_sample_name} \
                -o reg

            mkdir output_dir
            cp reg/*.vcf.gz* output_dir/

        else
            # Handle case where user provides PAR bed by running with --squash-ploidy over haploid region
            awk -v OFS="\t" '{print $1, 0, $2}' ~{reference.index} > genome_file.txt
            bedtools complement -i ~{par_bed} -g genome_file.txt -L > 'non-par.bed'
            bedtools complement -i non-par.bed -g genome_file.txt > 'regular-regions.bed'

            rtg format -o rtg_ref ~{reference.fasta}

            rtg vcfeval \
                ~{false="--all-records" true="" passing_only} \
                ~{true="--ref-overlap" false="" enable_ref_overlap} \
                --squash-ploidy \
                -b ~{base_vcf} \
                -c ~{query_vcf} \
                --bed-regions non-par.bed \
                ~{"-e " + evaluation_bed} \
                --vcf-score-field="~{score_field}" \
                --output-mode combine \
                --decompose \
                --roc-subset snp,indel \
                -t rtg_ref \
                ~{"--sample " + base_vcf_sample_name + "," + query_vcf_sample_name} \
                -o par

            rtg vcfeval \
                ~{false="--all-records" true="" passing_only} \
                ~{false="--squash-ploidy" true="" require_matching_genotypes} \
                ~{true="--ref-overlap" false="" enable_ref_overlap} \
                -b ~{base_vcf} \
                -c ~{query_vcf} \
                --bed-regions regular-regions.bed \
                ~{"-e " + evaluation_bed} \
                --vcf-score-field="~{score_field}" \
                --output-mode combine \
                --decompose \
                --roc-subset snp,indel \
                -t rtg_ref \
                ~{"--sample " + base_vcf_sample_name + "," + query_vcf_sample_name} \
                -o reg

            mkdir output_dir
            bcftools merge --force-samples "par/output.vcf.gz" "reg/output.vcf.gz" | bcftools sort -Oz -o "output_dir/output.vcf.gz"
            bcftools index -t "output_dir/output.vcf.gz"

        fi

        # Format ROC stats into table
        python3 << CODE
        import gzip
        import pandas as pd

        def parse_data(root_dir):
            full_df = pd.DataFrame()
            for Type in ['snp', 'indel']:
                file_path = f'{root_dir}/{Type}_roc.tsv.gz'

                header_lines = []
                # Read through file lines until hitting one without leading '#'
                with gzip.open(file_path, 'rt') as file:
                    for line in file:
                        if line[0] == '#':
                            header_lines += [line]
                        else:
                            break
                header_names = header_lines[-1].replace('#', '').replace('\n', '').split('\t')
                df = pd.read_csv(file_path, sep='\t', comment='#', header=None, names=header_names)

                rename_columns = {'score': 'Score', 'true_positives_baseline': 'TP_Base',
                      'false_positives': 'FP', 'true_positives_call': 'TP_Query', 'false_negatives': 'FN',
                      'precision': 'Precision', 'sensitivity': 'Recall', 'f_measure': 'F1_Score'}

                df = df.rename(columns=rename_columns)
                df['Type'] = Type.upper()
                df['Query_Name'] = "~{query_output_sample_name}"
                df['Base_Name'] = "~{base_output_sample_name}"

                full_df = pd.concat([full_df, df])

            return full_df

        reg_roc_summary = parse_data('reg')
        roc_summary = reg_roc_summary

        # If PAR bed file provided, also collect data from analysis over PAR region and combine stats
        if len("~{par_bed}") > 0:
            par_roc_summary = parse_data('par')
            merged_df = reg_roc_summary.merge(par_roc_summary, on=['Score', 'Type', 'Interval', 'Query_Name', 'Base_Name'], how='outer').fillna(0)
            for stat in ['TP_Base', 'FP', 'TP_Query', 'FN']:
                merged_df[stat] = merged_df[f'{stat}_x'] + merged_df[f'{stat}_y']
            merged_df['Precision'] = merged_df['TP_Query'] / (merged_df['TP_Query'] + merged_df['FP'])
            merged_df['Recall'] = merged_df['TP_Base'] / (merged_df['TP_Base'] + merged_df['FN'])
            merged_df['F1_Score'] = 2 * merged_df['Precision'] * merged_df['Recall'] / (merged_df['Precision'] + merged_df['Recall'])

            roc_summary = merged_df[
                ['Score', 'TP_Base', 'FP', 'TP_Query', 'FN', 'Precision', 'Recall', 'F1_Score', 'Type', 'Interval', 'Query_Name', 'Base_Name']
            ]

        roc_summary.to_csv('ROC_summary.tsv', sep='\t', index=False)

        CODE

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.0"
        preemptible: select_first([preemptible, 0])
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        cpu: runtimeAttributes.cpu
        memory: runtimeAttributes.memory + " GB"
    }

    output {
        File ROC_summary = "ROC_summary.tsv"

        File combined_output = "output_dir/output.vcf.gz"
        File combined_output_index = "output_dir/output.vcf.gz.tbi"
    }
}

task BCFToolsStats {
    input {
        File combined_vcfeval_output
        File combined_vcfeval_output_index

        File? stratifier_interval
        String stratifier_label

        String bcf_genotype
        String bcf_genotype_label

        String query_output_sample_name
        String base_output_sample_name

        String targets_overlap = "record"

        RuntimeAttributes runtimeAttributes = {"disk_size":2 * ceil(size(combined_vcfeval_output, "GB") + size(stratifier_interval, "GB")) + 10,
                                                  "cpu": 4, "memory": 8}
    }

    # Handle parsing bcf_genotype with surrounding quotes, with empty case handled separately
    String selection = if bcf_genotype!= "" then " && " + bcf_genotype + "'" else "'"

    String tp_base_selection = "-i 'INFO/BASE=\"TP\"" + selection
    String tp_query_selection = "-i 'INFO/CALL=\"TP\"" + selection
    String fp_selection = "-i '(INFO/CALL=\"FP\" || INFO/CALL=\"FP_CA\")" + selection
    String fn_selection = "-i '(INFO/BASE=\"FN\" || INFO/BASE=\"FN_CA\")" + selection
    String out_selection = "-i '(INFO/BASE=\"OUT\" || INFO/CALL=\"OUT\")" + selection
    String ign_selection = "-i '(INFO/BASE=\"IGN\" || INFO/CALL=\"IGN\")" + selection

    String targets_overlap_expression = " --targets-overlap " + targets_overlap

    command <<<
        set -xeuo pipefail

        # Subset combined output to each stat category; run in parallel
        # Note: subset to sample name and drop unused alts with -a BEFORE using above selector in stats command,
        # otherwise selecting for GT/TYPE all at once will apply to ALL sample/allele fields, giving misleading stats
        # Note: Use -T flag instead of -R for subsettng by region; the former streams the entire VCF file while latter does random
        # access lookup for each entry, which is extremely slow for even interval lists with just a few thousand entries
        # Also use --targets-overlap record to ensure events just overlapping the bed get counted (rather than being fully contained)
        bcftools view -s BASELINE --min-ac 1 -a -I ~{combined_vcfeval_output} | bcftools stats ~{tp_base_selection} ~{"-T " + stratifier_interval + targets_overlap_expression} > tp_base_stats.tsv &
        bcftools view -s CALLS --min-ac 1 -a -I ~{combined_vcfeval_output} | bcftools stats ~{tp_query_selection} ~{"-T " + stratifier_interval + targets_overlap_expression} > tp_query_stats.tsv &
        bcftools view -s CALLS --min-ac 1 -a -I ~{combined_vcfeval_output} | bcftools stats ~{fp_selection} ~{"-T " + stratifier_interval + targets_overlap_expression} > fp_stats.tsv &
        bcftools view -s BASELINE --min-ac 1 -a -I ~{combined_vcfeval_output} | bcftools stats ~{fn_selection} ~{"-T " + stratifier_interval + targets_overlap_expression} > fn_stats.tsv &
        bcftools stats ~{out_selection} -s- ~{"-T " + stratifier_interval + targets_overlap_expression} ~{combined_vcfeval_output} > out_stats.tsv &
        bcftools stats ~{ign_selection} -s- ~{"-T " + stratifier_interval + targets_overlap_expression} ~{combined_vcfeval_output} > ign_stats.tsv &

        # Wait for bcftools stats to finish for all selectors
        wait

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
            "TP_Query": make_dfs("tp_query_stats.tsv"),
            "FP": make_dfs("fp_stats.tsv"),
            "FN": make_dfs("fn_stats.tsv"),
            "OUT": make_dfs("out_stats.tsv"),
            "IGN": make_dfs("ign_stats.tsv")
        }

        # Make full SN (Summary Numbers) df
        full_sn_df = pd.DataFrame({'Category': df_dict['TP_Base']['SN_df']['Category']})
        for stat in df_dict:
            full_sn_df = df_dict[stat]['SN_df'].rename(columns={'Count': f'{stat}_Count'}).merge(full_sn_df, on='Category', how='outer')
        full_sn_df = full_sn_df.fillna(0)
        full_sn_df['Query_Name'] = "~{query_output_sample_name}"
        full_sn_df['Base_Name'] = "~{base_output_sample_name}"
        full_sn_df['Interval'] = "~{stratifier_label}"
        full_sn_df['BCF_Label'] = "~{bcf_genotype_label}"

        # Make full IDD (InDel Distribution) df
        full_idd_df = pd.DataFrame({'INDEL_Length': df_dict['TP_Base']['IDD_df']['INDEL_Length']})
        for stat in df_dict:
            full_idd_df = df_dict[stat]['IDD_df'].rename(columns={'Count': f'{stat}_Count'}).merge(full_idd_df, on='INDEL_Length', how='outer')
        full_idd_df = full_idd_df.fillna(0)
        full_idd_df['Query_Name'] = "~{query_output_sample_name}"
        full_idd_df['Base_Name'] = "~{base_output_sample_name}"
        full_idd_df['Interval'] = "~{stratifier_label}"
        full_idd_df['BCF_Label'] = "~{bcf_genotype_label}"

        # Make full ST (Substitution) df
        full_st_df = pd.DataFrame({'Substitution': df_dict['TP_Base']['ST_df']['Substitution']})
        for stat in df_dict:
            full_st_df = df_dict[stat]['ST_df'].rename(columns={'Count': f'{stat}_Count'}).merge(full_st_df, on='Substitution', how='outer')
        full_st_df = full_st_df.fillna(0)
        full_st_df['Query_Name'] = "~{query_output_sample_name}"
        full_st_df['Base_Name'] = "~{base_output_sample_name}"
        full_st_df['Interval'] = "~{stratifier_label}"
        full_st_df['BCF_Label'] = "~{bcf_genotype_label}"

        # Write files
        full_sn_df.to_csv('Full_SN.tsv', sep='\t', index=False)
        full_idd_df.to_csv('Full_IDD.tsv', sep='\t', index=False)
        full_st_df.to_csv('Full_ST.tsv', sep='\t', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.2"
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        cpu: runtimeAttributes.cpu
        memory: runtimeAttributes.memory + "GB"
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

        String? experiment
        String? extra_column_name
        String? extra_column_value

        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(size(ROC_summaries, "GB") + size(SN_summaries, "GB")
                                             + size(IDD_summaries, "GB") + size(ST_summaries, "GB")) + 10, "cpu": 2, "memory": 8}
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

        if "~{experiment}" != "":
            full_ROC['Experiment'] = "~{experiment}"

        if ("~{extra_column_name}" != "") and ("~{extra_column_value}" != ""):
            full_ROC["~{extra_column_name}"] = "~{extra_column_value}"

        full_ROC.to_csv('ROCStats.tsv', sep='\t', index=False)

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
            df['Precision'] = df['TP_Query'] / (df['TP_Query'] + df['FP'])
            df['Recall'] = df['TP_Base'] / (df['TP_Base'] + df['FN'])
            df['F1_Score'] = 2 * df['Precision'] * df['Recall'] / (df['Precision'] + df['Recall'])

        # Add optional labels
        for df in [full_SN, full_IDD, full_ST, simple_summary]:
            if "~{experiment}" != "":
                df['Experiment'] = "~{experiment}"

            if ("~{extra_column_name}" != "") and ("~{extra_column_value}" != ""):
                df["~{extra_column_name}"] = "~{extra_column_value}"

        # Reorder columns
        metadata_cols = ['Query_Name', 'Base_Name', 'Interval', 'Type']
        metadata_cols = metadata_cols + ["~{extra_column_name}"] if "~{extra_column_name}" != "" else metadata_cols
        metadata_cols = ['Experiment'] + metadata_cols if "~{experiment}" != "" else metadata_cols
        stat_cols = ['TP_Query', 'TP_Base', 'FP', 'FN', 'Precision', 'Recall', 'F1_Score', 'IGN', 'OUT']

        full_IDD = full_IDD[metadata_cols + ['INDEL_Type', 'INDEL_Length'] + stat_cols]
        full_ST = full_ST[metadata_cols + ['Substitution', 'Substitution_Type'] + stat_cols]
        simple_summary = simple_summary[metadata_cols + stat_cols]

        # Skip returning full_SN since it is redundant with simple_summary
        simple_summary.to_csv('SimpleSummary.tsv', sep='\t', index=False)
        full_IDD.to_csv('IndelDistributionStats.tsv', sep='\t', index=False)
        full_ST.to_csv('SNPSubstitutionStats.tsv', sep='\t', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        cpu: runtimeAttributes.cpu
        memory: runtimeAttributes.memory + "GB"
    }

    output {
        File simple_summary = "SimpleSummary.tsv"
        File IDD_combined_summaries = "IndelDistributionStats.tsv"
        File ST_combined_summaries = "SNPSubstitutionStats.tsv"
        File ROC_combined_summaries = "ROCStats.tsv"
    }
}
