version 1.0

struct RuntimeAttributes {
    Int disk_size
    Int cpu
    Int memory
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
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.1"
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        memory: runtimeAttributes.memory + " GB"
        cpu: runtimeAttributes.cpu
    }

    output {
        File output_table = "~{output_name}.tsv"
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
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.1"
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        memory: runtimeAttributes.memory + " GB"
        cpu: runtimeAttributes.cpu
    }

    output {
        File output_vcf = "output.vcf.gz"
        File output_vcf_index = "output.vcf.gz.tbi"
    }
}
