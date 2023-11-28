version 1.0

workflow CleanSVs {
    input {
        File input_vcf
        File? input_vcf_index

        String? output_name

        Int min_size = 50    # Minimum size of SV events

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
                min_size=min_size,
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
        File? input_vcf_index

        String output_name = "split"
        String output_suffix = "-splitMA"

        File? restriction_bed

        # Runtime parameters
        Int disk_size = ceil(2 * size(input_vcf, "GB")) + 100
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xe

        if [ -z "~{input_vcf_index}" ]; then
            bcftools index -t ~{input_vcf}
        fi

        bcftools norm -m -any -o "~{output_name}~{output_suffix}.vcf.gz" ~{"-R " + restriction_bed} ~{input_vcf}
        bcftools index -t -f "~{output_name}~{output_suffix}.vcf.gz"
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.0"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
    }

    output {
        File output_vcf = "~{output_name}~{output_suffix}.vcf.gz"
        File output_vcf_index = "~{output_name}~{output_suffix}.vcf.gz.tbi"
    }
}

# Used to add SVTYPE, SVLEN, and END fields to VCF containing just large INDELs
task AddAnnotations {
    input {
        File input_vcf
        File input_vcf_index

        String output_name = "annotated"
        String output_suffix = "-sv_annotated"

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
        bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%SVLEN\t%SVTYPE\n' truvari-SVs.vcf.gz > query.tsv

        # Process updated END coordinate from POS and SVLEN
        python3 << CODE
        import pandas as pd
        import numpy as np

        df = pd.read_csv('query.tsv', sep='\t', names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'SVLEN', 'SVTYPE'])
        df['SVLEN'] = df['SVLEN'].replace('.', 0).astype(int).apply(np.abs)
        # df['DIST'] = df.apply(lambda x: 1 if x['SVTYPE'] == 'INS' else x['SVLEN'], axis=1)
        # df['END'] = df['POS'] + df['DIST']
        df['END'] = df['POS'] + df['SVLEN']
        df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'END', 'SVLEN']].to_csv('annotations.tsv', sep='\t', index=False, header=False)
        CODE

        # Add annotations using bcftools
        # Start by removing old SVLEN from annotations; see this issue: https://github.com/samtools/bcftools/issues/1998
        bcftools annotate -x INFO/SVLEN truvari-SVs.vcf.gz -o truvari-SVs-noSVLEN.vcf.gz
        bgzip annotations.tsv
        tabix -s1 -b2 -e2 annotations.tsv.gz
        echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate in reference for SV">' > annotations.hdr
        echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length of ref and alt alleles">' >> annotations.hdr
        bcftools annotate -a annotations.tsv.gz -h annotations.hdr -c CHROM,POS,~ID,REF,ALT,INFO/END,INFO/SVLEN -o "~{output_name}~{output_suffix}.vcf.gz" truvari-SVs-noSVLEN.vcf.gz
        bcftools index -t -f "~{output_name}~{output_suffix}.vcf.gz"
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.0"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
    }

    output {
        File output_vcf = "~{output_name}~{output_suffix}.vcf.gz"
        File output_vcf_index = "~{output_name}~{output_suffix}.vcf.gz.tbi"
    }
}

# Replaces explicit INDEL sequences with abstractions, e.g. <DEL> or <INS>
# Also converts '.' FILTER to 'PASS' which is needed for Witty.er evaluation
task ConvertToAbstract {
    input {
        File input_vcf
        File input_vcf_index

        String output_name = "converted"
        String output_suffix = "-converted_abs"

        Int min_size = 50

        Boolean normalize_missing_filter_to_pass = true
        Boolean remove_bnd = true

        # Runtime parameters
        Int disk_size = ceil(2 * size(input_vcf, "GB")) + 100
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xe

        python3 << CODE
        import pysam
        import numpy as np

        with pysam.VariantFile("~{input_vcf}") as vcf:
            with pysam.VariantFile("~{output_name}~{output_suffix}.vcf.gz", "w", header=vcf.header) as output_vcf:
                for record in vcf:
                    ## Remove alt sequence to use abstract alleles; assumes split MA sites already
                    ref = record.alleles[0]
                    alt = record.alleles[1]

                    # Check relative difference in size meets size threshold or already abstract
                    # Also convert BND notation using ] or [ to abstract allele
                    if (np.abs(len(ref) - len(alt)) >= ~{min_size}) or (alt[0] == '<') or ('[' in alt) or (']' in alt):
                        # If SVTYPE has Number=A in header, record.info["SVTYPE"] would be a tuple
                        if isinstance(record.info['SVTYPE'], tuple):
                            svtype = record.info['SVTYPE'][0]    # Multiallelics already split, so just one entry
                        else:
                            svtype = record.info['SVTYPE']    # Should be string in this case

                        # Replace REF with its first char, and ALT with abstraction
                        # If SVTYPE is UNK, try to infer from alt alleles
                        if svtype == "UNK":
                            allele_type = alt.split('<')[-1].split('>')[0]
                            record.info["SVTYPE"] = allele_type
                        else:
                            new_alt = f'<{svtype}>' if alt[0] != '<' else alt
                            record.alleles = (record.alleles[0][0].upper(), new_alt)    # upper since Wittyer dislikes lowercase ref...

                        # Normalize FILTER to be PASS rather than missing
                        if (len(record.filter.values()) == 0) and ('~{normalize_missing_filter_to_pass}' == 'true'):
                            record.filter.clear()
                            record.filter.add('PASS')

                        # Check for BND and write only if not asked to remove
                        if not ( ('~{remove_bnd}' == 'true') and (svtype == 'BND') ):
                            if record.info['SVLEN'] > 0:    # sniffles actually will write variants with abstract alleles SVLEN = 0, which causes problem in Wittyer...
                                output_vcf.write(record)
        CODE

        bcftools index -t -f "~{output_name}~{output_suffix}.vcf.gz"
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/sv_docker:v1.0"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_ram + " GB"
        cpu: cpu
    }

    output {
        File output_vcf = "~{output_name}~{output_suffix}.vcf.gz"
        File output_vcf_index = "~{output_name}~{output_suffix}.vcf.gz.tbi"
    }
}
