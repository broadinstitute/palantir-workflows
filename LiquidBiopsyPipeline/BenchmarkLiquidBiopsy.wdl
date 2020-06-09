
workflow BenchmarkLiquidBiopsy {

   String bloodbiopsydocker
   String lego_plotter
   File reference
   File reference_index
   File reference_dict
   String reference_version

   String base_name

   File truth_vcf
   File truth_vcf_index
   File high_confidence_interval
   File high_confidence_bed
   File eval_vcf
   File eval_vcf_index

   File tumor_bam
   File tumor_bam_index
   File normal_bam
   File normal_bam_index
    
   String generate_plots_rscript

   call GenerateTruthStatistics {
      input:
         generate_plots_rscript = generate_plots_rscript,
         lego_plotter = lego_plotter,
         reference = reference,
         reference_index = reference_index,
         reference_dict = reference_dict,
         truth_vcf = truth_vcf,
         truth_vcf_index = truth_vcf_index,
         eval_vcf = eval_vcf, 
         eval_vcf_index = eval_vcf_index, 
         bam_file = tumor_bam,
         bam_index = tumor_bam_index,
         base_name = base_name,
         high_confidence_interval = high_confidence_interval,
         bloodbiopsydocker = bloodbiopsydocker
   }

   call CreatePlots {
      input: 
         bloodbiopsydocker = bloodbiopsydocker,
         base_name = base_name,
         truthTable = GenerateTruthStatistics.truth_table,
         evalTable = GenerateTruthStatistics.eval_table,
         pileup = GenerateTruthStatistics.pileup,
         generate_plots_rscript = generate_plots_rscript
   }

   call M2_Concordance {
      input: 
         bloodbiopsydocker = bloodbiopsydocker,
         truth_vcf = truth_vcf,
         truth_vcf_index = truth_vcf_index,
         eval_vcf = eval_vcf,
         eval_vcf_index = eval_vcf_index,
         intervals = high_confidence_interval, 
         base_name = base_name
   }

   Array[String] input_files = [truth_vcf,
                                eval_vcf,
                                GenerateTruthStatistics.TP_vcf,
                                GenerateTruthStatistics.FP_vcf,
                                GenerateTruthStatistics.FN_vcf, 
                                tumor_bam, 
                                normal_bam,
                                high_confidence_bed]
                                
   Array[String] input_names = ["Truth", "Eval Calls", "TP", "FP", "FN", "tumor", "normal", "high_confidence"]
   
   call GenerateIGVSession {
      input: 
         bloodbiopsydocker = bloodbiopsydocker,
         input_files = input_files,
         file_name =  base_name, 
         reference_version = reference_version, 
         input_names = input_names
   }
   output {
      File igv_session = GenerateIGVSession.igv_session
      File results_tsv = CreatePlots.results_tsv
      File summary_tsv = CreatePlots.summary_tsv
      File filter_analysis = M2_Concordance.filter_analysis
      File concordance_summary = M2_Concordance.summary
   }
}

# This task is used to assess the performance of the pipeline and will not be included in a final operational pipeline (since we won't know what the truth is).
task GenerateTruthStatistics {
   String generate_plots_rscript
   String bloodbiopsydocker
   String lego_plotter
   File reference
   File reference_index
   File reference_dict
   File truth_vcf
   File truth_vcf_index
   File high_confidence_interval
   File eval_vcf
   File eval_vcf_index
   File bam_file
   File bam_index
   String base_name

   Int disk_pad = 20
   Int disk_size = ceil(size(reference, "GB") +
                   size(reference_index, "GB") +
                   size(reference_dict, "GB") +
                   size(truth_vcf, "GB") +
                   size(truth_vcf_index, "GB") +
                   size(eval_vcf, "GB") +
                   size(eval_vcf_index, "GB") +
                   size(bam_file, "GB") +
                   size(bam_index, "GB") +
                   size(high_confidence_interval, "GB")) * 2 +
                   disk_pad

   Int? preemptible_attempts

   command <<<
      set -e

      # Calculate Pileup.  This is used later for R.
      gatk Pileup \
         -R ${reference} \
         -I ${bam_file} \
         -L ${truth_vcf}\
         -L ${eval_vcf} \
         -O ${base_name}.pileup

      # Subset variants to confident region and to just SNPs
      gatk SelectVariants -V ${eval_vcf} -L ${high_confidence_interval} -O confident.eval.vcf.gz -R ${reference} -select-type SNP
      gatk SelectVariants -V ${truth_vcf} -L ${high_confidence_interval} -O confident.truth.vcf.gz -R ${reference} -select-type SNP

      # Create VCFs of Truth Status (TP, FP, FN)
      gatk SelectVariants -V confident.eval.vcf.gz -O confident.eval.PASS.vcf.gz -R ${reference} --exclude-filtered
      gatk SelectVariants -V confident.truth.vcf.gz --concordance confident.eval.PASS.vcf.gz -O ${base_name}.TP.vcf.gz -R ${reference}
      gatk SelectVariants -V confident.eval.PASS.vcf.gz --discordance confident.truth.vcf.gz -O ${base_name}.FP.vcf.gz -R ${reference}
      gatk SelectVariants -V confident.truth.vcf.gz --discordance confident.eval.PASS.vcf.gz -O ${base_name}.FN.vcf.gz -R ${reference} -L confident.truth.vcf.gz

      # Create tables of Truth Status that will be read into R
      gatk VariantsToTable -F CHROM -F POS -F REF -F ALT -F FILTER -R ${reference} -O ${base_name}.truthTable.table -V confident.truth.vcf.gz
      gatk VariantsToTable -F CHROM -F POS -F REF -F ALT -F FILTER -F TLOD -R ${reference} -O ${base_name}.evalTable.table -V confident.eval.vcf.gz --show-filtered

      # Create lego plots
      # The FORMAT field has to be removed from the VCF in order for the legoplotter to work
      python2.7 ${lego_plotter} -m ${base_name}.TP.vcf.gz -t vcf -r ${reference} -b ${base_name}.TP.legoplot
      python2.7 ${lego_plotter} -m ${base_name}.FN.vcf.gz -t vcf -r ${reference} -b ${base_name}.FN.legoplot
      python2.7 ${lego_plotter} -m ${base_name}.FP.vcf.gz -t vcf -r ${reference} -b ${base_name}.FP.legoplot

   >>>
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: "6 GB"
      preemptible: select_first([preemptible_attempts, 0])
   }
   output {
      File truth_table = "${base_name}.truthTable.table"
      File eval_table = "${base_name}.evalTable.table"
      File pileup = "${base_name}.pileup"
      File TP_vcf = "${base_name}.TP.vcf.gz"
      File TP_vcf_idx = "${base_name}.TP.vcf.gz.tbi"
      File TP_legoplot = "${base_name}.TP.legoplot.pdf"
      File FP_vcf = "${base_name}.FP.vcf.gz"
      File FP_vcf_idx = "${base_name}.FP.vcf.gz.tbi"
      File FP_legoplot = "${base_name}.FP.legoplot.pdf"
      File FN_vcf = "${base_name}.FN.vcf.gz"
      File FN_vcf_idx = "${base_name}.FN.vcf.gz.tbi"
      File FN_legoplot = "${base_name}.FN.legoplot.pdf"
   }
}

task CreatePlots {
   String bloodbiopsydocker
   String generate_plots_rscript
   File truthTable
   File evalTable
   File pileup

   String base_name

   Int? preemptible_attempts

   Int disk_pad = 20
   Int disk_size = ceil(size(truthTable, "GB") +
                   size(evalTable, "GB") +
                   size(pileup, "GB")) * 2 +
                   disk_pad

   command <<<
      set -e
      Rscript -e "source('${generate_plots_rscript}');makeplot('${truthTable}', '${evalTable}', '${pileup}', '${base_name}')" 
   >>>

   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: "3 GB"
      preemptible: select_first([preemptible_attempts, 0])
   }
   output {
      File results_tsv = "${base_name}.results.tsv"
      File summary_tsv = "${base_name}.summary.results.tsv"

      File allele_fraction = "${base_name}.af.pdf"
      File binned_sensitivity = "${base_name}.binnedSensitivity.pdf"
      File call_status_vs_af = "${base_name}.callStatusVsAF.pdf"
      File depth = "${base_name}.depth.pdf"
      File depth_vs_alt = "${base_name}.depthVsAlt.pdf"
      File tlod = "${base_name}.tlod.pdf"
   }
}

task M2_Concordance {
   String bloodbiopsydocker
   String base_name
   File truth_vcf
   File truth_vcf_index
   File eval_vcf
   File eval_vcf_index
   File intervals

   Int? preemptible_attempts

   Int disk_pad = 20
   Int disk_size = ceil(size(truth_vcf, "GB") +
                   size(truth_vcf_index, "GB") +
                   size(eval_vcf, "GB") +
                   size(eval_vcf_index, "GB") +
                   size(intervals, "GB")) * 2 +
                   disk_pad

   command {

   gatk Concordance -L ${intervals} \
      -truth ${truth_vcf} \
      -eval ${eval_vcf}  \
      -tpfn tpfn.vcf \
      -ftnfn tnfn.vcf \
      -tpfp tpfp.vcf \
      -summary "${base_name}-summary.tsv" \
      --filter-analysis "${base_name}-filter-analysis.txt"
   }

   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: "3 GB"
      preemptible: select_first([preemptible_attempts, 0])
   }

   output {
      File filter_analysis = "${base_name}-filter-analysis.txt"
      File summary = "${base_name}-summary.tsv"
   }
}

# Creates an IGV session
# given a list of IGV compatible files (as strings).
# Reference is either "hg19" or "hg38".
task GenerateIGVSession {
   String bloodbiopsydocker
   Array[String] input_files
   String reference_version
   String file_name
   Array[String]? input_names
   Array[String] input_names_prefix = if defined(input_names) then prefix('-n ', select_first([input_names])) else []

   Int? preemptible_attempts
   Int disk_size = 20
    
   command {
      bash /usr/writeIGV.sh '${reference_version}' ${sep=" " input_files} ${sep=" " input_names_prefix}  > "${file_name}.xml"
   }
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: "3 GB"
      preemptible: select_first([preemptible_attempts, 0])
   }
   output {
      File igv_session = "${file_name}.xml"
   }
}



