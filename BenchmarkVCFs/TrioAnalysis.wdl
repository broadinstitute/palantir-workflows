version 1.0

workflow TrioAnalysis {
    input {
        File child_vcf
        File child_vcf_index
        String sex_of_child

        File father_vcf
        File father_vcf_index

        File mother_vcf
        File mother_vcf_index

        File ref_fasta
        File ref_index

        Array[File]? bed_files
        Array[String]? bed_labels
    }

    call ComputeMendelianConcordance {
        input:
            child_vcf = child_vcf,
            child_vcf_index = child_vcf_index,
            sex_of_child = sex_of_child,
            father_vcf = father_vcf,
            father_vcf_index = father_vcf_index,
            mother_vcf = mother_vcf,
            mother_vcf_index = mother_vcf_index,
            ref_fasta = ref_fasta,
            ref_index = ref_index
    }

    if (defined(bed_files) && defined(bed_labels)) {
        call AnnotateVcfRegions {
            input:
                input_vcf = ComputeMendelianConcordance.mendelian_output,
                input_vcf_index = ComputeMendelianConcordance.mendelian_output_index,
                bed_files = select_first([bed_files]),
                bed_labels = select_first([bed_labels])
        }
    }

    call MakeSummaryTable {
        input:
            annotated_vcf = select_first([AnnotateVcfRegions.annotated_vcf, ComputeMendelianConcordance.mendelian_output]),
            annotated_vcf_index = select_first([AnnotateVcfRegions.annotated_vcf_index, ComputeMendelianConcordance.mendelian_output_index]),
            bed_labels = bed_labels
    }

    output {
        File mendelian_output_vcf = select_first([AnnotateVcfRegions.annotated_vcf, ComputeMendelianConcordance.mendelian_output])
        File mendelian_output_vcf_index = select_first([AnnotateVcfRegions.annotated_vcf_index, ComputeMendelianConcordance.mendelian_output_index])

        File mendelian_violation_table = MakeSummaryTable.mendelian_violation_table
        File mendelian_uncertainty_table = MakeSummaryTable.mendelian_uncertainty_table
        File mendelian_violation_counts = MakeSummaryTable.mendelian_violation_counts
        File mendelian_uncertainty_counts = MakeSummaryTable.mendelian_uncertainty_counts
        File mendelian_violation_cdf = MakeSummaryTable.mendelian_violation_cdf
        File mendelian_uncertainty_cdf = MakeSummaryTable.mendelian_uncertainty_cdf
    }

}

# Following task computes the mendelian concordance stats using rtg mendelian
task ComputeMendelianConcordance {
    input {
        File child_vcf
        File child_vcf_index
        String sex_of_child

        File father_vcf
        File father_vcf_index

        File mother_vcf
        File mother_vcf_index

        File ref_fasta
        File ref_index

        Boolean lenient = false
    }

    String lenient_command = if (lenient) then "--lenient" else ""

    command <<<
        set -xueo pipefail

        # Merge into trio
        bcftools merge -m all -o trio-merge.vcf.gz ~{child_vcf} ~{father_vcf} ~{mother_vcf}

        # Re-normalized/left align (sometimes collapsing in merge can shift things)
        bcftools norm -m +any -o trio-norm.vcf.gz -f ~{ref_fasta} trio-merge.vcf.gz

        # Add pedigree info to header for rtg
        CHILD_NAME=$(bcftools query -l ~{child_vcf})
        FATHER_NAME=$(bcftools query -l ~{father_vcf})
        MOTHER_NAME=$(bcftools query -l ~{mother_vcf})

        echo "##PEDIGREE=<Child=${CHILD_NAME}, Father=${FATHER_NAME}, Mother=${MOTHER_NAME}>" > add_header.txt
        echo "##SAMPLE=<ID=${CHILD_NAME},Sex=~{sex_of_child}>" >> add_header.txt
        bcftools annotate -h add_header.txt trio-norm.vcf.gz -o trio.vcf.gz -Wtbi

        # Make rtg reference
        rtg format -o rtg_ref ~{ref_fasta}

        # Run rtg mendelian
        rtg mendelian \
            -i trio.vcf.gz \
            -t rtg_ref \
            -o mendelian_output.vcf.gz \
            ~{lenient_command}
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
        disks: "local-disk " + 200 + " HDD"
        cpu: 8
        memory: 16 + " GB"
    }

    output {
        File mendelian_output = "mendelian_output.vcf.gz"
        File mendelian_output_index = "mendelian_output.vcf.gz.tbi"
    }
}

# Use bedtools and bcftools to annotate a VCF using array of BED files with labels
task AnnotateVcfRegions {
    input {
        File input_vcf
        File input_vcf_index
        Array[File] bed_files
        Array[String] bed_labels

        Int disk_size = ceil(2 * size(input_vcf, "GB")) + 100
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xeuo pipefail

        bed_files=( ~{sep=" " bed_files} )
        bed_labels=( ~{sep=" " bed_labels} )

        # Move input VCF to avoid overwriting in loop
        mv ~{input_vcf} input.vcf.gz
        mv ~{input_vcf_index} input.vcf.gz.tbi

        # Annotate VCF with each BED file and corresponding label
        for i in ${!bed_files[@]}; do
            bed_file=${bed_files[$i]}
            bed_label=${bed_labels[$i]}

            bgzip $bed_file
            tabix -p bed "${bed_file}.gz"
            bcftools annotate input.vcf.gz \
                -a "${bed_file}.gz" \
                -c CHROM,POS \
                -h <(echo "##INFO=<ID=${bed_label},Number=0,Type=Flag,Description=\"Contained in ${bed_label}\">") \
                -m "+${bed_label}" \
                -o annotated.vcf.gz
            mv annotated.vcf.gz input.vcf.gz
            bcftools index -t input.vcf.gz
        done

        # Rename final output
        mv input.vcf.gz annotated_output.vcf.gz
        mv input.vcf.gz.tbi annotated_output.vcf.gz.tbi
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
    }

    output {
        File annotated_vcf = "annotated_output.vcf.gz"
        File annotated_vcf_index = "annotated_output.vcf.gz.tbi"
    }
}

# A task for making summary table of Mendelian violation counts
task MakeSummaryTable {
    input {
        File annotated_vcf
        File annotated_vcf_index

        Array[String]? bed_labels

        String format_value = "GQ"
        Int default_format_value = 0
        String format_agg = "min"  # One of: min, max, mean, median, sum

        String? experiment
    }

    String initial_query_tab = if (defined(bed_labels)) then "\\t%" else ""

    command <<<
        set -xueo pipefail

        bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%MCV\t%MCU[\t%GT][\t%~{format_value}]~{initial_query_tab}~{sep="\\t%" bed_labels}\n' ~{annotated_vcf} > mendelian_table-preheader.tsv
        echo -e "CHROM\tPOS\tREF\tALT\tMCV\tMCU\tChild_GT\tFather_GT\tMother_GT\tChild_~{format_value}\tFather_~{format_value}\tMother_~{format_value}~{true="\\t" false="" defined(bed_labels)}~{sep="\\t" bed_labels}" > mendelian_table_header.tsv

        cat mendelian_table_header.tsv mendelian_table-preheader.tsv > mendelian_table.tsv

        python3 << CODE
        import pandas as pd
        import numpy as np

        intervals = ["~{sep="\", \"" bed_labels}"]
        df = pd.read_csv("mendelian_table.tsv", sep="\t", na_values=["."])

        # Fill in missing values with default
        format_cols = [c for c in df.columns if '_~{format_value}' in c]
        df[format_cols] = df[format_cols].fillna(~{default_format_value})
        df[intervals] = df[intervals].fillna(0)

        mv_df = df[~df['MCV'].isna()]
        mcu_df = df[~df['MCU'].isna()]

        mv_counts_df = mv_df.value_counts(intervals).reset_index(name='MCV')
        mcu_counts_df = mcu_df.value_counts(intervals).reset_index(name='MCU')

        def make_cdf(df):
            df['Agg_~{format_value}'] = df.apply(lambda row: np.~{format_agg}([row['Child_~{format_value}'], row['Father_~{format_value}'], row['Mother_~{format_value}']]), axis=1)

            agg_counts = df['Agg_~{format_value}'].value_counts().reset_index().sort_values(by='Agg_~{format_value}', ascending=False)
            agg_counts['~{format_agg}_~{format_value}_cdf'] = agg_counts['count'].cumsum()
            return agg_counts

        mv_cdf = make_cdf(mv_df)
        mv_cdf['Interval'] = 'Whole Genome'
        for i in intervals:
            new_cdf = make_cdf(mv_df[mv_df[i] == 1])
            new_cdf['Interval'] = i
            mv_cdf = pd.concat([mv_cdf, new_cdf])

        mcu_cdf = make_cdf(mcu_df)
        mcu_cdf['Interval'] = 'Whole Genome'
        for i in intervals:
            new_cdf = make_cdf(mcu_df[mcu_df[i] == 1])
            new_cdf['Interval'] = i
            mcu_cdf = pd.concat([mcu_cdf, new_cdf])

        if '~{experiment}':
            for d in [mv_df, mcu_df, mv_counts_df, mcu_counts_df, mv_cdf, mcu_cdf]:
                d['Experiment'] = "~{experiment}"

        mv_df.to_csv("mendelian_violation_table.tsv.gz", sep="\t", index=False)
        mcu_df.to_csv("mendelian_uncertainty_table.tsv.gz", sep="\t", index=False)
        mv_counts_df.to_csv("mendelian_violation_counts.tsv.gz", sep="\t", index=False)
        mcu_counts_df.to_csv("mendelian_uncertainty_counts.tsv.gz", sep="\t", index=False)
        mv_cdf.to_csv("mendelian_violation_cdf.tsv.gz", sep="\t", index=False)
        mcu_cdf.to_csv("mendelian_uncertainty_cdf.tsv.gz", sep="\t", index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
        disks: "local-disk " + 200 + " HDD"
        cpu: 2
        memory: 4 + " GB"
    }

    output {
        File mendelian_violation_table = "mendelian_violation_table.tsv.gz"
        File mendelian_uncertainty_table = "mendelian_uncertainty_table.tsv.gz"
        File mendelian_violation_counts = "mendelian_violation_counts.tsv.gz"
        File mendelian_uncertainty_counts = "mendelian_uncertainty_counts.tsv.gz"
        File mendelian_violation_cdf = "mendelian_violation_cdf.tsv.gz"
        File mendelian_uncertainty_cdf = "mendelian_uncertainty_cdf.tsv.gz"
    }
}