version 1.0

task DetectHPVIntegrationBreakpoints {
    input {
        String output_basename
        File bam
        File bai

        Int cpu = 1
        Int memory_gb = 16
        Int disk_size_gb = ceil((2.5 * size(bam, "GiB")) + 50)
    }

    command <<<
        python /breakpoint_detector-v3.7.py ~{output_basename} ~{bam} .
    >>>

    output {
        File analysis_log = "~{output_basename}_analysis_log.txt"
        File breakpoints = "~{output_basename}_breakpoints.txt"
        File detailed_integration_summary = "~{output_basename}_detailed_integration_summary.txt"
        File integration_breakpoints = "~{output_basename}_integration_breakpoints.txt"
        File integration_summary = "~{output_basename}_integration_summary.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/breakpoint_detector@sha256:6795ccbcfa825c39cdd5a11ddac89673c445ef73266fcd7a7b0fb0f7f0001550"
    }
}

# Check HPV16_Ref coverage
# Make MSA & phyml tree
task Sublineages {
    input {
        String output_basename
        File bam
        File bai
        File reference
        File hpv16_sublineages
        Int read_threshold = 10

        Int cpu = 4
        Int memory_gb = 32
        Int disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
    }

    command <<<
        READS=$(samtools view -c ~{bam} HPV16_Ref || echo 0)

        if [ "$READS" -lt ~{read_threshold} ]; then
            LINE="~{output_basename},~{output_basename},HPV16_negative,NA"
            echo "run_id,library_id,closest_sublineage,patristic_distance" > ~{output_basename}.sublineage_call.csv
            echo "$LINE" >> ~{output_basename}.sublineage_call.csv
            echo "Only $READS reads (<~{read_threshold}); marking HPV16_negative."
            touch ~{output_basename}.combo.afa ~{output_basename}.combo.phy ~{output_basename}.combo.phy_phyml_stats.txt ~{output_basename}.combo.phy_phyml_tree.txt
        else
            bcftools mpileup --max-depth 8000 -Ou -f ~{reference} -r HPV16_Ref ~{bam} | \
            bcftools call -Ou -mv | \
            bcftools norm -f ~{reference} -Oz -o ~{output_basename}.vcf.gz

            tabix ~{output_basename}.vcf.gz

            samtools faidx ~{reference} HPV16_Ref | \
            bcftools consensus ~{output_basename}.vcf.gz -o ~{output_basename}.consensus.fasta

            cat ~{output_basename}.consensus.fasta ~{hpv16_sublineages} > ~{output_basename}.combo.fasta
            muscle -threads $(nproc) -align ~{output_basename}.combo.fasta -output ~{output_basename}.combo.afa

            seqret -osformat2 phylip -sequence ~{output_basename}.combo.afa -outseq ~{output_basename}.combo.phy

            phyml -i ~{output_basename}.combo.phy

            touch ~{output_basename}.sublineage_call.csv
        fi
    >>>

    output {
        File multiple_sequence_alignment = "~{output_basename}.combo.afa"
        File phylip_formatted_msa = "~{output_basename}.combo.phy"
        File phylogenetic_tree_stats = "~{output_basename}.combo.phy_phyml_stats.txt"
        File phylogenetic_tree = "~{output_basename}.combo.phy_phyml_tree.txt"
        File sublineage_call = "~{output_basename}.sublineage_call.csv"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hpv_sublineage@sha256:10f9ff9e349c440dd616b3b94f9e6928b0d4e822c78cabbe8aacfaf08e804e5f"
    }
}

# Extract nearest neighbor via Toytree and write to sample CSV
# Draw tree with Toytree
task SublineagesDrawTree {
    input {
        String output_basename
        File phylogenetic_tree
        File sublineage_call_in

        Int cpu = 1
        Int memory_gb = 32
        Int disk_size_gb = 64
    }

    command <<<
if [ -s "~{phylogenetic_tree}" ]; then
python <<CODE
import toytree
import toyplot.pdf

t = toytree.tree("~{phylogenetic_tree}")
df = t.distance.get_tip_distance_matrix(df = True)
d = df.loc["HPV16_Ref"].drop("HPV16_Ref")

with open("~{output_basename}.sublineage_call.csv", 'w') as f:
    f.write("run_id,library_id,closest_sublineage,patristic_distance\n")
    f.write("~{output_basename}" + "," + "~{output_basename}" + "," + str(d.idxmin()) + "," + "{:.8f}".format(d.min()))

    canvas = toytree.tree("~{phylogenetic_tree}").draw(node_labels = False)[0]
    toyplot.pdf.render(canvas, "~{output_basename}.combo.phy_phyml_tree.pdf")
CODE
else
    touch ~{output_basename}.combo.phy_phyml_tree.pdf
    cat ~{sublineage_call_in} > ~{output_basename}.sublineage_call.csv
fi
    >>>

    output {
        File phylogenetic_tree_visualization = "~{output_basename}.combo.phy_phyml_tree.pdf"
        File sublineage_call = "~{output_basename}.sublineage_call.csv"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hpv_toytree@sha256:9b056932890befab08c7c96e4cec884f205ceca1c7aac819d66204ea96854b7f"
    }
}

task HPVHighRiskSNPs {
    input {
        String output_basename
        File bam
        File bai
        File high_risk_snps_hpv
        File reference
        String genotype = "HPV16_Ref"

        Int cpu = 2
        Int memory_gb = 16
        Int disk_size_gb = ceil((2.5 * size(bam, "GiB")) + 50)
    }

    command <<<
        bcftools mpileup -Ou --fasta-ref ~{reference} --max-depth 800000 --regions ~{genotype} --annotate FORMAT/AD,INFO/AD,FORMAT/AD,FORMAT/ADF ~{bam} \
        | bcftools call -Ou -mv -o ~{output_basename}.vcf.gz

        echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tDP\tINFO/AC\tINFO/AD\tAF" | cat - <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%DP\t%INFO/AC\t%INFO/AD\n' ~{output_basename}.vcf.gz \
        | awk -F'\t' '{split($8,AD,","); if (AD[1]+AD[2] == 0) AF=0; else AF=AD[2]/(AD[1]+AD[2]); print $0 "\t" AF }') > ~{output_basename}.variants.txt

        awk -F'\t' 'FNR==NR{a[$1]=$1;next} {if($2 in a) {b[FILENAME]=$2}} END{for(i in b) {print i"\t"b[i]}}' ~{high_risk_snps_hpv} ~{output_basename}.variants.txt > ~{output_basename}.hr_snps_found.txt
    >>>

    output {
        File high_risk_snps_found = "~{output_basename}.hr_snps_found.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.4"
    }
}

workflow HPVDeepSeekTertiaryAnalysis {
    input {
        String output_basename
        File tumor_bam
        File tumor_bai
        File high_risk_snps_hpv
        File reference
        File reference_fai
        File reference_dict
        File hpv16_sublineages
    }

    call DetectHPVIntegrationBreakpoints {
        input:
            bam = tumor_bam,
            bai = tumor_bai,
            output_basename = output_basename
    }

    call Sublineages {
        input:
            bam = tumor_bam,
            bai = tumor_bai,
            reference = reference,
            hpv16_sublineages = hpv16_sublineages,
            output_basename = output_basename
    }

    call SublineagesDrawTree {
        input:
            phylogenetic_tree = Sublineages.phylogenetic_tree,
            sublineage_call_in = Sublineages.sublineage_call,
            output_basename = output_basename
    }

    call HPVHighRiskSNPs {
        input:
            bam = tumor_bam,
            bai = tumor_bai,
            high_risk_snps_hpv = high_risk_snps_hpv,
            reference = reference,
            output_basename = output_basename
    }

    output {
        File analysis_log = DetectHPVIntegrationBreakpoints.analysis_log
        File breakpoints = DetectHPVIntegrationBreakpoints.breakpoints
        File detailed_integration_summary = DetectHPVIntegrationBreakpoints.detailed_integration_summary
        File integration_breakpoints = DetectHPVIntegrationBreakpoints.integration_breakpoints
        File integration_summary = DetectHPVIntegrationBreakpoints.integration_summary
        File multiple_sequence_alignment = Sublineages.multiple_sequence_alignment
        File phylip_formatted_msa = Sublineages.phylip_formatted_msa
        File phylogenetic_tree_stats = Sublineages.phylogenetic_tree_stats
        File phylogenetic_tree = Sublineages.phylogenetic_tree
        File phylogenetic_tree_visualization = SublineagesDrawTree.phylogenetic_tree_visualization
        File sublineage_call = SublineagesDrawTree.sublineage_call
        File high_risk_snps_found = HPVHighRiskSNPs.high_risk_snps_found
    }
}