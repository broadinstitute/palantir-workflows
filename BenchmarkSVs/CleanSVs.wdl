version 1.0

workflow CleanSVs {
    input {
        File input_vcf
        File input_vcf_index

        String? output_name

        Boolean split_MA = true
        Boolean add_annotations = true
        Boolean convert_to_abstract = true
    }

    if (split_MA) {
        call SplitMASites {
            input:
                input_vcf=input_vcf,
                input_vcf_index=input_vcf_index,
                output_name=output_name
        }
    }

    if (add_annotations) {
        call AddAnnotations {
            input:
                input_vcf=select_first([SplitMASites.output_vcf, input_vcf]),
                input_vcf_index=select_first([SplitMASites.output_vcf_index, input_vcf_index]),
                output_name=output_name
        }
    }

    if (convert_to_abstract) {
        call ConvertToAbstract {
            input:
                input_vcf=select_first([AddAnnotations.output_vcf, SplitMASites.output_vcf, input_vcf]),
                input_vcf_index=select_first([AddAnnotations.output_vcf_index, SplitMASites.output_vcf_index, input_vcf_index]),
                output_name=output_name
        }
    }

    output {
        File output_vcf = select_first([ConvertToAbstract.output_vcf, AddAnnotations.output_vcf, SplitMASites.output_vcf, input_vcf])
        File output_vcf_index = select_first([ConvertToAbstract.output_vcf_index, AddAnnotations.output_vcf_index, SplitMASites.output_vcf_index, input_vcf_index])
    }
}

# Splits MA sites into biallelics
task SplitMASites {
    input {
        File input_vcf
        File input_vcf_index

        String output_name = "split"

        # Runtime parameters
        Int disk_size = ceil(2 * size(input_vcf, "GB")) + 100
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xe

        bcftools norm -m -any -o "~{output_name}.vcf.gz" ~{input_vcf}
        bcftools index -t -f "~{output_name}.vcf.gz"
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.5"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
    }

    output {
        File output_vcf = "~{output_name}.vcf.gz"
        File output_vcf_index = "~{output_name}.vcf.gz.tbi"
    }
}

# Used to add SVTYPE, SVLEN, and END fields to VCF containing just large INDELs
task AddAnnotations {
    input {
        File input_vcf
        File input_vcf_index

        String output_name = "annotated"

        Int min_size = 50

        # Runtime parameters
        Int disk_size = ceil(2 * size(input_vcf, "GB")) + 100
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xe

        # Add SVTYPE and SVLEN fields
        truvari anno svinfo -m ~{min_size} -o truvari-annotated.vcf.gz ~{input_vcf}
        bcftools view -i 'INFO/SVTYPE!="."' -o truvari-SVs.vcf.gz truvari-annotated.vcf.gz
        bcftools index -t -f truvari-SVs.vcf.gz

        # Add END for INDELs using bcftools because pysam doesn't handle these well
        # See: https://github.com/pysam-developers/pysam/issues/1200
        bcftools query -f'%CHROM\t%POS\t%POS\t%SVLEN\t%SVTYPE\n' truvari-SVs.vcf.gz > query.tsv

        # Process updated END coordinate from POS and SVLEN
        python3 << CODE
        import pandas as pd
        import numpy as np

        df = pd.read_csv('query.tsv', sep='\t', names=['CHROM', 'POS', 'POS2', 'SVLEN', 'SVTYPE'])
        df['DIST'] = df.apply(lambda x: 1 if x['SVTYPE'] == 'INS' else np.abs(x['SVLEN']), axis=1)
        df['END'] = df['POS'] + df['DIST']
        df[['CHROM', 'POS', 'POS2', 'END']].to_csv('annotations.tsv', sep='\t', index=False, header=False)
        CODE

        # Add annotations using bcftools
        bgzip annotations.tsv
        tabix -s1 -b2 -e3 annotations.tsv.gz
        echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate in reference for SV">' > annotations.hdr
        bcftools annotate -a annotations.tsv.gz -h annotations.hdr -c CHROM,FROM,TO,INFO/END -o "~{output_name}.vcf.gz" truvari-SVs.vcf.gz
        bcftools index -t -f "~{output_name}.vcf.gz"
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.5"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
    }

    output {
        File output_vcf = "~{output_name}.vcf.gz"
        File output_vcf_index = "~{output_name}.vcf.gz.tbi"
    }
}

# Replaces explicit INDEL sequences with abstractions, e.g. <DEL> or <INS>
# Also converts '.' FILTER to 'PASS' which is needed for Witty.er evaluation
task ConvertToAbstract {
    input {
        File input_vcf
        File input_vcf_index

        String output_name = "converted"

        # Runtime parameters
        Int disk_size = ceil(2 * size(input_vcf, "GB")) + 100
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xe

        python3 << CODE
        import pysam

        with pysam.VariantFile("~{input_vcf}") as vcf:
            with pysam.VariantFile("~{output_name}.vcf.gz", "w", header=vcf.header) as output_vcf:
                for record in vcf:
                    # Remove alt sequence to use abstract alleles; assumes split MA sites already
                    # Replace REF with its first char, and ALT with abstraction
                    record.alleles = (record.alleles[0][0], f'<{record.info["SVTYPE"]}>')
                    if len(record.filter.values()) == 0:
                        record.filter.clear()
                        record.filter.add('PASS')
                    output_vcf.write(record)
        CODE

        bcftools index -t -f "~{output_name}.vcf.gz"
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.5"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
    }

    output {
        File output_vcf = "~{output_name}.vcf.gz"
        File output_vcf_index = "~{output_name}.vcf.gz.tbi"
    }
}
