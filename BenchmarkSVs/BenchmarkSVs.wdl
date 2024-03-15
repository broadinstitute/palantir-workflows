version 1.0

import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/Utilities/WDLs/CreateIGVSession.wdl" as IGV

struct RuntimeAttributes {
    Int disk_size
    Int cpu
    Int memory
}

workflow BenchmarkSVs {
    input {
        File? base_vcf
        File? base_vcf_index
        Array[String] base_sample_names

        File comp_vcf
        File comp_vcf_index
        Array[String] comp_sample_names

        String? experiment

        File ref_fasta
        File ref_fai

        Boolean perform_qc = true
        Boolean run_truvari = true

        File? evaluation_bed
        Float? evaluation_pct    # Defaults to checking at least one base overlap evaluation_bed

        Array[File] bed_regions = []
        Array[String] bed_labels = []
        Int breakpoint_padding = 20

        Array[Int] svlen_bin_cutoffs = [100, 250, 1000, 2500, 10000, 25000]

        Boolean create_igv_session = true
        Array[File]? optional_igv_bams
    }

    # Subset to evaluation regions for remainder of analysis
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

    # Replace input files with subset versions
    File subset_comp_vcf = select_first([SubsetComp.output_vcf, comp_vcf])
    File subset_comp_vcf_index = select_first([SubsetComp.output_vcf_index, comp_vcf_index])
    File subset_base_vcf = select_first([SubsetBase.output_vcf, base_vcf])
    File subset_base_vcf_index = select_first([SubsetBase.output_vcf_index, base_vcf_index])

    # QC Tasks
    if (perform_qc) {
        String qc_comp_output_name = if (defined(base_vcf)) then "qc_stats-comp" else "combined_qc_stats"
        call CollectQcMetrics as CompQcMetrics {
            input:
                input_vcf=subset_comp_vcf,
                input_vcf_index=subset_comp_vcf_index,
                sample_name=comp_sample_names[0],
                experiment=experiment,
                label_suffix="Comp",
                ref_fai=ref_fai,
                bed_regions=bed_regions,
                bed_labels=bed_labels,
                svlen_bin_cutoffs=svlen_bin_cutoffs
        }

        if (length(bed_regions) > 0) {
            call AddIntervalOverlapStats as CompQcMetricsIntervals {
                input:
                    input_table=CompQcMetrics.stats,
                    table_label="qc_comp",
                    output_name=qc_comp_output_name,
                    bed_regions=bed_regions,
                    bed_labels=bed_labels,
                    ref_fai=ref_fai,
                    breakpoint_padding=breakpoint_padding
            }
        }

        if (defined(base_vcf)) {
            call CollectQcMetrics as BaseQcMetrics {
                input:
                    input_vcf=subset_base_vcf,
                    input_vcf_index=subset_base_vcf_index,
                    sample_name=base_sample_names[0],
                    experiment=experiment,
                    label_suffix="Base",
                    ref_fai=ref_fai,
                    bed_regions=bed_regions,
                    bed_labels=bed_labels,
                    svlen_bin_cutoffs=svlen_bin_cutoffs
            }

            if (length(bed_regions) > 0) {
                call AddIntervalOverlapStats as BaseQcMetricsIntervals {
                    input:
                        input_table=BaseQcMetrics.stats,
                        table_label="qc_base",
                        output_name="qc_base",
                        bed_regions=bed_regions,
                        bed_labels=bed_labels,
                        ref_fai=ref_fai,
                        breakpoint_padding=breakpoint_padding
                }
            }

            call CombineSummaries as CombineQcMetrics {
                input:
                    tables=[select_first([CompQcMetricsIntervals.output_table, CompQcMetrics.stats]),
                            select_first([BaseQcMetricsIntervals.output_table, BaseQcMetrics.stats])],
                    output_file_name="combined_qc_stats.tsv"
            }
        }

        File qc_file = select_first([CombineQcMetrics.combined_file, CompQcMetrics.stats])
    }

    # Benchmarking - Truvari tasks
    ## Run Truvari over each sample name pair
    if (defined(base_vcf) && run_truvari) {
        scatter(sample_pair in zip(base_sample_names, comp_sample_names)) {
            call RunTruvari {
                input:
                    base_vcf=subset_base_vcf,
                    base_vcf_index=subset_base_vcf_index,
                    base_sample_name=sample_pair.left,
                    comp_vcf=subset_comp_vcf,
                    comp_vcf_index=subset_comp_vcf_index,
                    comp_sample_name=sample_pair.right,
                    experiment=experiment
            }

            call AddIntervalOverlapStats as TruvariTpBaseIntervals {
                input:
                    input_table=RunTruvari.tp_base_table,
                    table_label="truvari_bench",
                    output_name=sample_pair.left+"tp_base_intervals",
                    bed_regions=bed_regions,
                    bed_labels=bed_labels,
                    ref_fai=ref_fai,
                    breakpoint_padding=breakpoint_padding
            }

            call AddIntervalOverlapStats as TruvariTpCompIntervals {
                input:
                    input_table=RunTruvari.tp_comp_table,
                    table_label="truvari_bench",
                    output_name=sample_pair.right+"tp_comp_intervals",
                    bed_regions=bed_regions,
                    bed_labels=bed_labels,
                    ref_fai=ref_fai,
                    breakpoint_padding=breakpoint_padding
            }

            call AddIntervalOverlapStats as TruvariFpIntervals {
                input:
                    input_table=RunTruvari.fp_table,
                    table_label="truvari_bench",
                    output_name=sample_pair.right+"fp_intervals",
                    bed_regions=bed_regions,
                    bed_labels=bed_labels,
                    ref_fai=ref_fai,
                    breakpoint_padding=breakpoint_padding
            }

            call AddIntervalOverlapStats as TruvariFnIntervals {
                input:
                    input_table=RunTruvari.fn_table,
                    table_label="truvari_bench",
                    output_name=sample_pair.left+"fn_intervals",
                    bed_regions=bed_regions,
                    bed_labels=bed_labels,
                    ref_fai=ref_fai,
                    breakpoint_padding=breakpoint_padding
            }

            call ComputeTruvariIntervalSummaryStats {
                input:
                    base_sample_name=sample_pair.left,
                    comp_sample_name=sample_pair.right,
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
                    base_sample_name=sample_pair.left,
                    comp_sample_name=sample_pair.right,
                    experiment=experiment,
                    ref_fai=ref_fai
            }
        }

        call CombineSummaries as CombineTruvariSummaries {
            input:
                tables=ComputeTruvariIntervalSummaryStats.full_summary,
                output_file_name="truvari_bench_summary.tsv"
        }

        call CombineSummaries as CombineFPClosest {
            input:
                tables=CollectTruvariClosestStats.fp_closest,
                output_file_name="truvari_fp_closest.tsv"
        }

        call CombineSummaries as CombineFNClosest {
            input:
                tables=CollectTruvariClosestStats.fn_closest,
                output_file_name="truvari_fn_closest.tsv"
        }

        if (create_igv_session) {
            call IGV.CreateIGVSession as IGVSession {
                input:
                    bams=optional_igv_bams,
                    vcfs=flatten([select_all(RunTruvari.tp_base), select_all(RunTruvari.tp_comp),
                                    select_all(RunTruvari.fp), select_all(RunTruvari.fn)]),
                    interval_files=bed_regions,
                    reference=ref_fasta
            }
        }
    }

    Array[File] qc_files = select_all([qc_file])
    call CombineFiles {
        input:
            qc_files=qc_files,
            truvari_files=select_all([CombineTruvariSummaries.combined_file, CombineFPClosest.combined_file, CombineFNClosest.combined_file])
    }

    output {
        # QC outputs
        File? qc_summary = qc_file

        # Truvari outputs
        File? truvari_bench_summary = CombineTruvariSummaries.combined_file
        File? truvari_fp_closest = CombineFPClosest.combined_file
        File? truvari_fn_closest = CombineFNClosest.combined_file

        # Compressed/combined outputs
        File? combined_files = CombineFiles.combined_files

        # IGV Session
         File? igv_session = IGVSession.igv_session
    }
}

task SubsetEvaluation {
    input {
        File input_vcf
        File input_vcf_index

        File evaluation_bed
        Float? evaluation_pct

        # Runtime parameters
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(input_vcf, "GB")) + 100, "cpu": 4, "memory": 16}
    }

    command <<<
        set -xueo pipefail

        # Subset input VCF to only sites which intersect evaluation regions over pct overlap threshold
        bcftools view -h ~{input_vcf} > header.txt
        bedtools intersect -a ~{input_vcf} -b ~{evaluation_bed} ~{"-f " + evaluation_pct} -u > variants.vcf

        cat header.txt variants.vcf | bcftools view -o output.vcf.gz -
        bcftools index -t output.vcf.gz

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.0"
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        memory: runtimeAttributes.memory + " GB"
        cpu: runtimeAttributes.cpu
    }

    output {
        File output_vcf = "output.vcf.gz"
        File output_vcf_index = "output.vcf.gz.tbi"
    }
}

task AddIntervalOverlapStats {
    input {
        File input_table
        String table_label

        String output_name

        Array[File] bed_regions
        Array[String] bed_labels

        File ref_fai

        Int breakpoint_padding = 20

        # Runtime Parameters
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(input_table, "GB")) + 100, "cpu": 4, "memory": 16}
    }

    command <<<
        set -xeuo pipefail

        # Clean ref_fai into genome file for bedtools
        cut -f1,2 ~{ref_fai} > ref.genome

        # Collect interval bed overlap stats
        # Assumes input_table will have first three columns: CHROM, POS0, END0 with no header
        generate_interval_stats () {
            INTERVAL_LABEL=$1
            INTERVAL_FILE=$2
            INPUT_LABEL=$3
            VCF_FILE=$4

            echo -e "${INTERVAL_LABEL}-count\t${INTERVAL_LABEL}-overlap" > header.txt
            # WARNING: annotate does not preserve order of input bed regions!
            cut -f1-3 $VCF_FILE > vcf.bed
            bedtools annotate -both -i vcf.bed -files $INTERVAL_FILE | bedtools sort -i - | rev | cut -f1-2 | rev > "${INTERVAL_LABEL}-${INPUT_LABEL}-preheader.bed"
            cat header.txt "${INTERVAL_LABEL}-${INPUT_LABEL}-preheader.bed" > "${INTERVAL_LABEL}-${INPUT_LABEL}-annotated.bed"

            # Also add in breakpoint overlap stats using POS and END values
            echo -e "${INTERVAL_LABEL}-LBEND-count\t${INTERVAL_LABEL}-LBEND-overlap" > pos-header.txt
            echo -e "${INTERVAL_LABEL}-RBEND-count\t${INTERVAL_LABEL}-RBEND-overlap" > end-header.txt
            awk '{OFS="\t"; print $1, $2, $2}' vcf.bed | bedtools slop -b ~{breakpoint_padding} -i - -g ref.genome > pos.bed
            awk '{OFS="\t"; print $1, $3, $3}' vcf.bed | bedtools slop -b ~{breakpoint_padding} -i - -g ref.genome > end.bed
            bedtools annotate -both -i pos.bed -files $INTERVAL_FILE | bedtools sort -i - | rev | cut -f1-2 | rev > "${INTERVAL_LABEL}-${INPUT_LABEL}-pos-preheader.bed"
            bedtools annotate -both -i end.bed -files $INTERVAL_FILE | bedtools sort -i - | rev | cut -f1-2 | rev > "${INTERVAL_LABEL}-${INPUT_LABEL}-end-preheader.bed"
            cat pos-header.txt "${INTERVAL_LABEL}-${INPUT_LABEL}-pos-preheader.bed" > "${INTERVAL_LABEL}-${INPUT_LABEL}-pos-annotated.bed"
            cat end-header.txt "${INTERVAL_LABEL}-${INPUT_LABEL}-end-preheader.bed" > "${INTERVAL_LABEL}-${INPUT_LABEL}-end-annotated.bed"

            # Paste together three overlap stats files
            paste "${INTERVAL_LABEL}-${INPUT_LABEL}-annotated.bed" "${INTERVAL_LABEL}-${INPUT_LABEL}-pos-annotated.bed" "${INTERVAL_LABEL}-${INPUT_LABEL}-end-annotated.bed" > "${INTERVAL_LABEL}-${INPUT_LABEL}-full-annotated.bed"
        }

        INTERVAL_FILES=(~{sep=' ' bed_regions})
        INTERVAL_LABELS=(~{sep=' ' bed_labels})

        # split header from table content
        # head -1 "~{input_table}" > table.hdr
        tail -n+2 "~{input_table}" > table_content.tsv

        for i in {0..~{length(bed_regions)-1}}
        do
            INTERVAL_LABEL="${INTERVAL_LABELS[$i]}"
            CURRENT_FILE="${INTERVAL_FILES[$i]}"
            generate_interval_stats $INTERVAL_LABEL $CURRENT_FILE "~{table_label}" table_content.tsv
        done

        # Combine across interval beds
        paste "~{input_table}" *-~{table_label}-full-annotated.bed > "~{output_name}.tsv"

        # Correct POS0/END0 back to original VCF POS/END coordinates
        python3 << CODE
        import pandas as pd

        df = pd.read_csv('~{output_name}.tsv', sep='\t')
        df['POS'] = df['POS'] + 1
        df['END'] = df['END'] + 1
        df.to_csv('~{output_name}.tsv', sep='\t', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.0"
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        memory: runtimeAttributes.memory + " GB"
        cpu: runtimeAttributes.cpu
    }

    output {
        File output_table = "~{output_name}.tsv"
    }
}

task CollectQcMetrics {
    input {
        File input_vcf
        File input_vcf_index

        String sample_name = ""
        String experiment = ""
        String label_suffix = "Comp"

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
            tagged.vcf.gz | awk -v OFS='\t' '{ print $0, "~{sample_name}-~{label_suffix}", "~{experiment}-~{label_suffix}" }' >> qc_stats-pre_intervals.tsv

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

task RunTruvari {
    input {
        File base_vcf
        File base_vcf_index
        String? base_sample_name

        File comp_vcf
        File comp_vcf_index
        String? comp_sample_name

        String experiment = "Experiment"

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

        Boolean debug_mode = false

        # Runtime parameters
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(base_vcf, "GB") + 2 * size(comp_vcf, "GB")) + 100, "cpu": 8, "memory": 16}
    }

    command <<<
        set -xueo pipefail

        # Subset to sites for given sample names since Truvari flags don't seem to work right
        if [ -z "~{base_sample_name}~{comp_sample_name}" ]
        then
            mv ~{base_vcf} base.vcf.gz
            mv ~{base_vcf_index} base.vcf.gz.tbi
            mv ~{comp_vcf} comp.vcf.gz
            mv ~{comp_vcf_index} comp.vcf.gz.tbi
        else
            bcftools view ~{"-s " + base_sample_name} --min-ac 1 -o base.vcf.gz ~{base_vcf}
            bcftools index -t base.vcf.gz
            bcftools view ~{"-s " + comp_sample_name} --min-ac 1 -o comp.vcf.gz ~{comp_vcf}
            bcftools index -t comp.vcf.gz
        fi

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

        # Use Python to collect the output results and compile across intervals into one table
        python3 << CODE
        import pandas as pd
        import json

        with open(f'output_dir/summary.json') as file:
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

        bcftools query -i 'INFO/SVTYPE!="."' -f"${MAIN_QUERY}" output_dir/tp-base.vcf.gz | bedtools sort -i - > tp-base-preheader.bed
        cat truvari_header.txt tp-base-preheader.bed > tp-base.bed

        bcftools query -i 'INFO/SVTYPE!="."' -f"${MAIN_QUERY}" output_dir/fn.vcf.gz | bedtools sort -i - > fn-preheader.bed
        cat truvari_header.txt fn-preheader.bed > fn.bed

        bcftools query -i 'INFO/SVTYPE!="."' -f"${MAIN_QUERY}" output_dir/tp-comp.vcf.gz | bedtools sort -i - > tp-comp-preheader.bed
        cat truvari_header.txt tp-comp-preheader.bed > tp-comp.bed

        bcftools query -i 'INFO/SVTYPE!="."' -f"${MAIN_QUERY}" output_dir/fp.vcf.gz | bedtools sort -i - > fp-preheader.bed
        cat truvari_header.txt fp-preheader.bed > fp.bed

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.0"
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

task CombineSummaries {
    input {
        Array[File] tables
        String output_file_name
    }

    Int disk_size = ceil(2 * size(tables, "GB")) + 50

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
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.0"
        disks: "local-disk " + disk_size + " HDD"
        memory: 8 + " GB"
        cpu: 4
    }

    output {
        File combined_file = "~{output_file_name}"
    }
}

task CombineFiles {
    input {
        Array[File] qc_files
        Array[File] truvari_files
    }

    command <<<
        set -xueo pipefail

        mkdir benchmark_outputs

        mkdir benchmark_outputs/qc_files
        if [ ~{length(qc_files)} -gt 0 ]
        then
            cp ~{sep=" " qc_files} benchmark_outputs/qc_files/
        fi

        mkdir benchmark_outputs/truvari_files
        if [ ~{length(truvari_files)} -gt 0 ]
        then
            cp ~{sep=" " truvari_files} benchmark_outputs/truvari_files/
        fi

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