version 1.0

workflow ConcordanceHail {
    input {
        Array[File] vcf_files
        Array[String] truth_sample_names
        String basename
        File truth_vcf
        File sites_csv
    }

    scatter (i in range(length(vcf_files))) {
        call ComputeConcordance {
            input:
                vcf = vcf_files[i],
                truth_sample_name = truth_sample_names[i],
                truth_vcf = truth_vcf,
                sites_csv = sites_csv
        }
    }

    call ConcatDataFrames {
        input:
            data_frames = ComputeConcordance.concordance_report,
            basename = basename
    }

    output {
        File concordance_report = ConcatDataFrames.concatenated_data_frame
    }
}

task ComputeConcordance {
    input {
        File vcf
        String truth_sample_name
        File truth_vcf
        File sites_csv
    }

    parameter_meta {
        vcf: {
            localization_optional : true
        }
        truth_vcf: {
            localization_optional : true
        }
    }

    command <<<
        set -euo pipefail
        cat << "EOF" > concordance.py
import hail as hl
import numpy as np
import pandas as pd

hl.init(default_reference='GRCh38')
truth_vcf = hl.import_vcf("~{truth_vcf}", force_bgz=True)
rg37 = hl.get_reference('GRCh37')
rg38 = hl.get_reference('GRCh38')  
rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)
sites_pd = pd.read_csv("~{sites_csv}")

sites_ht = hl.Table.from_pandas(sites_pd[['chromosome','position_b37']].sort_values(['chromosome','position_b37']).query('position_b37!=-1'))
lifted = hl.liftover(hl.locus(sites_ht.chromosome, sites_ht.position_b37,'GRCh37'),'GRCh38')
sites_ht = sites_ht.annotate(contig_hg38 = lifted.contig,
                             position_hg38 = lifted.position,
                            locus = lifted)
sites_ht = sites_ht.key_by(sites_ht.locus).cache()

truth_vcf_gsa_sites = truth_vcf.filter_rows(hl.is_defined(sites_ht[truth_vcf.locus]))
imputed_vcf = hl.import_vcf("~{vcf}",
                           force_bgz=True)
imputed_vcf_gsa_sites = imputed_vcf.filter_rows(hl.is_defined(sites_ht[imputed_vcf.locus]))
imputed_vcf_gsa_sites = imputed_vcf_gsa_sites.annotate_cols(wgs_name="~{truth_sample_name}")
imputed_vcf_gsa_sites = imputed_vcf_gsa_sites.key_cols_by(imputed_vcf_gsa_sites.wgs_name)
_, sample_conc, _ = hl.concordance(imputed_vcf_gsa_sites, truth_vcf_gsa_sites)
conc_np = np.array(sample_conc.concordance.collect()[0])
concordance = (sum([conc_np[i][i] for i in range(5)]) + conc_np[2][0])/(conc_np.sum()-conc_np[0].sum() - conc_np[2][1] - 
                                                                        conc_np[3][1] - conc_np[4][1])
concorance_df = pd.DataFrame({'sample_name':["~{truth_sample_name}"],
                                 'concordance':[concordance]})
concorance_df.to_csv("~{truth_sample_name}_concordance.tsv", sep="\t", index=False)
EOF

        python3 concordance.py
    >>>

    output {
        File concordance_report = "~{truth_sample_name}_concordance.tsv"
    }

    runtime {
        docker: "hailgenetics/hail:0.2.126-py3.11"
        disks: "local-disk 100 HDD"
        memory: "32 GiB"
        cpu: 8   
        preemptible: 0 
    }
}

task ConcatDataFrames {
    input {
        Array[File] data_frames
        String basename
    }

    command <<<
        python3 << "EOF"
        import pandas as pd
        dfs = [pd.read_csv(f, sep="\t") for f in ["~{sep='","' data_frames}"]]
        concat_df = pd.concat(dfs)
        concat_df.to_csv("~{basename}_concordance.tsv", sep="\t", index=False)
        EOF
    >>>

    output {
        File concatenated_data_frame = "~{basename}_concordance.tsv"
    }

    runtime {
        docker : "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk 100 HDD"
    }
}