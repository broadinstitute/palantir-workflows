version 1.0

workflow Glimpse2ImputationQC {
    input {
        File imputed_vcf
        
        Int preemptible = 1
        String docker = "hailgenetics/hail:0.2.126-py3.11"
        Int cpu = 4
        Int mem_gb = 16
        File? monitoring_script
    }

    call Glimpse2ImputationQCTask {
        input:
            imputed_vcf = imputed_vcf,
            preemptible = preemptible,
            docker = docker,
            cpu = cpu,
            mem_gb = mem_gb,
            monitoring_script = monitoring_script
    }

    output {
        File qc_metrics = Glimpse2ImputationQCTask.qc_metrics
    }
}

task Glimpse2ImputationQCTask {
    input {
        File imputed_vcf
        
        Int preemptible = 1
        String docker = "hailgenetics/hail:0.2.126-py3.11"
        Int cpu = 4
        Int mem_gb = 16
        File? monitoring_script
    }

    parameter_meta {
        imputed_vcf: {
                        localization_optional: true
                    }
    }

    Int disk_size_gb = 100
    
    command <<<
        set -euo pipefail

        ~{"bash " + monitoring_script + " > monitoring.log &"}

        cat <<'EOF' > script.py
import hail as hl
hl.init(default_reference='GRCh38', idempotent=True)
vcf = hl.import_vcf('~{imputed_vcf}', force_bgz=True)
qc = hl.sample_qc(vcf)
qc.cols().export('qc.tsv')
EOF
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File qc_metrics = "qc.tsv"
        File? monitoring = "monitoring.log"
    }
}
