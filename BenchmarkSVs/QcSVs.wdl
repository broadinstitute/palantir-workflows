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
        File vcf
        File vcf_index

        String? experiment

        File ref_fasta
        File ref_fai

        File? evaluation_bed
        Float? evaluation_pct    # Defaults to checking at least one base overlap evaluation_bed

        Array[File] bed_regions = []
        Array[String] bed_labels = []
        Int breakpoint_padding = 20

        Array[Int] svlen_bin_cutoffs = [100, 250, 1000, 2500, 10000, 25000]
    }

    # Subset to evaluation regions for remainder of analysis
    if (defined(evaluation_bed)) {
        call SV_tasks.SubsetEvaluation as SubsetComp {
            input:
                input_vcf=vcf,
                input_vcf_index=vcf_index,
                evaluation_bed=select_first([evaluation_bed]),
                evaluation_pct=evaluation_pct
        }
    }

    # Replace input files with subset versions
    File subset_vcf = select_first([SubsetComp.output_vcf, vcf])
    File subset_vcf_index = select_first([SubsetComp.output_vcf_index, vcf_index])

    # QC Tasks
    call CollectQcMetrics as QcMetrics {
        input:
            input_vcf=subset_vcf,
            input_vcf_index=subset_vcf_index,
            experiment=experiment,
            ref_fai=ref_fai,
            bed_regions=bed_regions,
            bed_labels=bed_labels,
            svlen_bin_cutoffs=svlen_bin_cutoffs
    }

    if (length(bed_regions) > 0) {
        call SV_tasks.AddIntervalOverlapStats as QcMetricsIntervals {
            input:
                input_table=QcMetrics.stats,
                table_label="qc_stats",
                output_name="qc_stats",
                bed_regions=bed_regions,
                bed_labels=bed_labels,
                ref_fai=ref_fai,
                breakpoint_padding=breakpoint_padding
        }
    }

    File qc_file = select_first([QcMetricsIntervals.output_table, QcMetrics.stats])

    output {
        File qc_summary = qc_file
    }
}

task CollectQcMetrics {
    input {
        File input_vcf
        File input_vcf_index

        String experiment = ""

        File ref_fai

        Array[File] bed_regions = []
        Array[String] bed_labels = []
        Int breakpoint_padding = 20

        Array[Int] svlen_bin_cutoffs = [100, 250, 1000, 2500, 10000, 25000]
        Array[Int] af_bin_cutoffs = [1, 10, 50]
        Boolean count_singletons_separate = true

        # Runtime parameters
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(input_vcf, "GB")) + 100, "cpu": 4, "memory": 16}
    }

    command <<<
        set -xueo pipefail

        # Add allele count annotations, in case they're missing
        bcftools +fill-tags -o tagged.vcf.gz "~{input_vcf}"
        bcftools index -t -f tagged.vcf.gz

        # Main query formatting
        # Use POS0/END0 for coherence with bed coordinates later
        MAIN_HEADER="CHROM\tPOS\tEND\tQUAL\tFILTER\tSVTYPE\tSVLEN\tAC\tAC_Het\tAC_Hom\tAN\tExcHet\tHWE\tAF\tMAF\tNS\tSample\tExperiment"
        MAIN_QUERY="%CHROM\t%POS0\t%END0\t%QUAL\t%FILTER\t%INFO/SVTYPE\t%INFO/SVLEN\t%AC\t%AC_Het\t%AC_Hom\t%AN\t%ExcHet\t%HWE\t%AF\t%MAF\t%NS\n"

        # Create table of relevant stats from variant annotations
        # Requires input to have INFO fields: END, SVTYPE, SVLEN
        echo -e "${MAIN_HEADER}" > qc_stats-pre_intervals.tsv
        bcftools query -i 'INFO/SVTYPE!="."' \
            -f"${MAIN_QUERY}" \
            tagged.vcf.gz | awk -v OFS='\t' '{ print $0, "~{experiment}" }' >> qc_stats-pre_intervals.tsv

        # Add extra stats to QC
        python3 << CODE
        import pandas as pd

        qc_df = pd.read_csv('qc_stats-pre_intervals.tsv', sep='\t')

        # Add AC_Ref counts
        qc_df['AC_Ref'] = qc_df['NS'] - qc_df['AC_Het'] - qc_df['AC_Hom'] / 2

        # Bin AF
        AF_Bin_cutoffs = [~{default="" sep=", " af_bin_cutoffs}]
        AF_Bin_start = [f'<{AF_Bin_cutoffs[0]}%']
        AF_Bins_middle = [f'{AF_Bin_cutoffs[i]}-{AF_Bin_cutoffs[i+1]}%' for i in range(len(AF_Bin_cutoffs)-1)]
        AF_Bin_end = [f'>{AF_Bin_cutoffs[-1]}%']
        AF_Bins = AF_Bin_start + AF_Bins_middle + AF_Bin_end

        def bin_af(af):
            for i, pct in enumerate(AF_Bin_cutoffs):
                if af < pct/100:
                    return AF_Bins[i]
            return AF_Bins[-1]

        # Fill missing AF with 0
        qc_df['AF_Bin'] = qc_df['AF'].replace('.', 0).astype(float).apply(bin_af)
        if "~{count_singletons_separate}" == "true":
            qc_df['AF_Bin'] = qc_df.apply(lambda x: 'AC=1' if x['AC'] == 1 else x['AF_Bin'], axis=1)

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

        qc_df['SVLEN_Bin'] = qc_df['SVLEN'].replace('.', 0).apply(bin_svlen)

        qc_df.to_csv('qc_stats-pre_intervals.tsv', sep='\t', index=False)

        CODE

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.0"
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        memory: runtimeAttributes.memory + " GB"
        cpu: runtimeAttributes.cpu
    }

    output {
        File stats = "qc_stats-pre_intervals.tsv"
    }
}
