
workflow GenerateDuplexConsensusBams {

   # reference files 
   File reference 
   File reference_index 
   File reference_bwa_index 
   File reference_pac 
   File reference_amb 
   File reference_ann 
   File reference_bwt 
   File reference_sa 
   File reference_dict 
   File reference_alt
   File known_indels 
   File known_indels_index 
   File variant_eval_gold_standard 
   File variant_eval_gold_standard_index  
   File dbsnp 
   File dbsnp_index
   String bwa_path 

   # sample specific files 
   File bam_file 
   File bam_index  
   String base_name 
   File target_intervals 
   File bait_intervals 

   # parameters 
   Int minimum_base_quality 
   Int allowable_umi_distance 
   String minimum_consensus_reads
   String min_reads
   Float frac_Ns
   String? consensus_extra_filter_args
   Float? downsample_probability
   Int num_clip_bases_five_prime
   Int? num_clip_bases_three_prime
   Boolean? run_bwa_mem_on_raw
   Boolean run_bwa_mem_on_raw_or_default = select_first([run_bwa_mem_on_raw, false])
   Int compression_level

   # scripts
   String process_duplex_coverage_rscript
   String calcFS_script

   # runtime 
   Int? preemptible_attempts
   ## Use as a last resort to increase the disk given to every task in case of ill behaving data
   Int? emergency_extra_disk
   String bloodbiopsydocker
   
   # This is added to every task as padding, should increase if systematically you need more disk for every call
   Int disk_pad = 10 + select_first([emergency_extra_disk,0])

   # Get the version of BWA that we are using.
   call GetBwaVersion {
      input:
         bloodbiopsydocker=bloodbiopsydocker,
         bwa_path = bwa_path, 
         preemptible_attempts = preemptible_attempts
   }

   if (defined(downsample_probability)) {
      call DownsampleSam {
         input:
            bloodbiopsydocker = bloodbiopsydocker, 
            bam_file = bam_file,
            bam_index = bam_index, 
            downsample_probability = downsample_probability,
            base_name = base_name, 
            preemptible_attempts = preemptible_attempts,
            disk_pad = disk_pad
      }
   }

   if (run_bwa_mem_on_raw_or_default) {

      call QuerySortSam {
         input:
            bloodbiopsydocker = bloodbiopsydocker, 
            input_bam = select_first([DownsampleSam.output_bam, bam_file]),
            base_name = base_name,
            preemptible_attempts = preemptible_attempts,
            disk_pad = disk_pad
      }

      call AlignRawBamWithBwaMem  {
         input: 
            input_bam = QuerySortSam.bam,
            bwa_version = GetBwaVersion.version,
            output_bam_basename = base_name,
            ref_fasta = reference,
            ref_fasta_index = reference_index,
            ref_dict = reference_dict,
            ref_alt = reference_alt,
            ref_amb = reference_amb,
            ref_ann = reference_ann,
            ref_bwt = reference_bwt,
            ref_pac = reference_pac,
            ref_sa = reference_sa,
            compression_level = compression_level,
            preemptible_attempts = preemptible_attempts, 
            disk_pad = disk_pad
      }
   }

   File preprocessed_raw_bam = select_first([AlignRawBamWithBwaMem.output_bam, DownsampleSam.output_bam, bam_file])
   File preprocessed_raw_bam_index = select_first([AlignRawBamWithBwaMem.output_bam_index, DownsampleSam.output_bam_index, bam_index])

   # Collect HS or Targeted PCR metrics after deduplication by start and stop
   # position (but not incluing UMIs).
   call CollectSelectionMetrics as BeforeUmiDeduplicationSelectionMetrics {
      input:
         bloodbiopsydocker=bloodbiopsydocker,
         bam_file = preprocessed_raw_bam,
         bam_index = preprocessed_raw_bam_index,
         base_name = base_name,
         reference = reference,
         reference_index = reference_index,
         reference_dict = reference_dict,
         target_intervals = target_intervals,
         bait_intervals = bait_intervals,
         suffix = ".noumi",
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   }

   # Assigns groups of reads in duplicate sets
   # into potentially finer grained duplicate sets
   # using their UMIs
   call FGBioGroupReadsByUmi {
      input:
         bloodbiopsydocker=bloodbiopsydocker,
         bam_file = preprocessed_raw_bam,
         bam_index = preprocessed_raw_bam_index,
         base_name = base_name,
         allowable_umi_distance = allowable_umi_distance,
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   }

   # Collect metrics on duplicate sets
   call CollectDuplexSeqMetrics {
      input:
         bloodbiopsydocker = bloodbiopsydocker,
         target_intervals = target_intervals,
         bam_file = FGBioGroupReadsByUmi.output_bam,
         base_name = base_name,
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   }

   call CalculateDuplexMetrics {
      input: 
         bloodbiopsydocker  = bloodbiopsydocker,
         duplex_family_sizes = CollectDuplexSeqMetrics.duplex_family_sizes, 
         calcFS_script = calcFS_script,
         preemptible_attempts = preemptible_attempts
   }

   # Calls duplex consensus reads generated from the same 
   # double-stranded source molecule. The min reads argument
   # refers to the minimum number of input reads to a consensus 
   # read. The first values applies to the final consensus read, 
   # the second value to one single-strand consensus, and the 
   # the last value to the other single-strand consensus. 
   call CallDuplexConsensusReads {
      input:
         bloodbiopsydocker=bloodbiopsydocker,
         bam_file = FGBioGroupReadsByUmi.output_bam,
         base_name = base_name,
         min_reads = min_reads, 
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   }

   call AlignReadsAndMBA as AlignReadsAndMBAConsensus {
      input:
         bloodbiopsydocker=bloodbiopsydocker,
         bam_file = CallDuplexConsensusReads.output_bam,
         reference = reference,
         reference_bwa_index = reference_bwa_index,
         reference_pac = reference_pac,
         reference_amb = reference_amb,
         reference_ann = reference_ann,
         reference_bwt = reference_bwt,
         reference_sa = reference_sa,
         reference_dict = reference_dict,
         bwa_version = GetBwaVersion.version,
         bwa_path = bwa_path,
         base_name = base_name,
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   }

   call FilterConsensusReads as FilterDuplexConsensusReads {
      input:
         bloodbiopsydocker=bloodbiopsydocker,
         bam_file = AlignReadsAndMBAConsensus.output_bam,
         bam_index = AlignReadsAndMBAConsensus.output_bam_index,
         reference = reference,
         base_name = base_name,
         minimum_consensus_reads = minimum_consensus_reads,
         minimum_base_quality = minimum_base_quality, 
         frac_Ns = frac_Ns,
         consensus_extra_filter_args = consensus_extra_filter_args, 
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   }

   call ClipBam {
      input:
         bloodbiopsydocker = bloodbiopsydocker,
         bam_file = FilterDuplexConsensusReads.output_bam,
         bam_index =  FilterDuplexConsensusReads.output_bam_index,
         base_name = base_name,
         num_clip_bases_five_prime =  num_clip_bases_five_prime,
         num_clip_bases_three_prime = num_clip_bases_three_prime,
         reference = reference,
         reference_index = reference_index, 
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   } 


   # Recompute base quality scores after consensus calls
   # are made.
   call BQSRWithoutBinning as BQSRDuplex {
      input:
         bloodbiopsydocker=bloodbiopsydocker,
         bam_file = ClipBam.output_bam,
         bam_index = ClipBam.output_bam_index,
         reference = reference,
         reference_index = reference_index,
         reference_dict = reference_dict,
         dbsnp = dbsnp,
         dbsnp_index= dbsnp_index,
         known_indels = known_indels,
         known_indels_index = known_indels_index,
         variant_eval_gold_standard = variant_eval_gold_standard,
         variant_eval_gold_standard_index = variant_eval_gold_standard_index,
         base_name = base_name,
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   }

   # Summarize base quality score statistics using BQSR after consensus
   # calling.
   call QualityScoreDistribution as QualityScoreDistributionDuplexConsensus  {
      input:
         bloodbiopsydocker=bloodbiopsydocker,
         bam_file = BQSRDuplex.output_bam,
         bam_index = BQSRDuplex.output_bam_index,
         base_name = base_name,
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   }

   # Collect HS or Targeted PCR metrics after conesnsus calls
   # for the tumor and normal samples.
   call CollectSelectionMetrics as DuplexSelectionMetrics {
      input:
         bloodbiopsydocker=bloodbiopsydocker,
         bam_file = BQSRDuplex.output_bam,
         bam_index = BQSRDuplex.output_bam_index,
         base_name = base_name,
         reference = reference,
         reference_index = reference_index,
         reference_dict = reference_dict,
         target_intervals = target_intervals,
         bait_intervals = bait_intervals,
         suffix = ".duplex", 
         preemptible_attempts = preemptible_attempts, 
         disk_pad = disk_pad
   }

   # Collect duplex fragment depth of coverage.
   call CollectDepthOfCoverage as CollectDuplexDepthOfCoverage {
      input: 
         interval_list = target_intervals,
         reference = reference,
         reference_index = reference_index,
         reference_dict = reference_dict,
         bam_file = BQSRDuplex.output_bam,
         bam_index = BQSRDuplex.output_bam_index,
         base_name = base_name,
         extra_arguments = "--countType COUNT_FRAGMENTS_REQUIRE_SAME_BASE",
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   }
   
   # Collect raw fragment depth of coverage.
   # This task can take as much as 20 hours, so it is recommended to set
   # preemptible_attempts to 0.
   call CollectDepthOfCoverage as CollectRawStartStopDepthOfCoverage {
      input: 
         interval_list = target_intervals,
         reference = reference,
         reference_index = reference_index,
         reference_dict = reference_dict,
         bam_file = preprocessed_raw_bam,
         bam_index = preprocessed_raw_bam_index,
         base_name = base_name,
         extra_arguments = "--countType COUNT_FRAGMENTS",
         preemptible_attempts = 0,
         disk_pad = disk_pad
   }
   
   # Collect raw read depth of coverage.
   # This task can take as much as 20 hours, so it is recommended to set
   # preemptible_attempts to 0.
   call CollectDepthOfCoverage as CollectRawReadDepthOfCoverage {
      input: 
         interval_list = target_intervals,
         reference = reference,
         reference_index = reference_index,
         reference_dict = reference_dict,
         bam_file = preprocessed_raw_bam,
         bam_index = preprocessed_raw_bam_index,
         base_name = base_name,
         extra_arguments = "--countType COUNT_READS -drf DuplicateRead",
         preemptible_attempts = 0,
         disk_pad = disk_pad
   }

   call CollectStatisticsByCoverage {
      input: 
         bloodbiopsydocker = bloodbiopsydocker,
         process_duplex_coverage_rscript = process_duplex_coverage_rscript,
         base_name = base_name,
         raw_depth = CollectRawReadDepthOfCoverage.depth_of_coverage,
         start_stop_depth = CollectRawStartStopDepthOfCoverage.depth_of_coverage,
         duplex_depth = CollectDuplexDepthOfCoverage.depth_of_coverage,
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   } 
   
   output {
      
      File duplex_output_bam = BQSRDuplex.output_bam
      File duplex_output_bam_index = BQSRDuplex.output_bam_index
      File raw_output_bam = preprocessed_raw_bam
      File raw_output_bam_index = preprocessed_raw_bam_index

      File duplex_selection_metrics = DuplexSelectionMetrics.output_selection_metrics
      File duplex_per_target_selection_metrics  = DuplexSelectionMetrics.output_per_target_selection_metrics 
      File duplex_theoretical_sensitivity  = DuplexSelectionMetrics.output_theoretical_sensitivity 

      Int mean_duplex_depth = CollectStatisticsByCoverage.mean_duplex_depth
      Int mean_raw_depth = CollectStatisticsByCoverage.mean_raw_depth
      
      File duplex_family_sizes = CollectDuplexSeqMetrics.duplex_family_sizes

      Float mean_family_size = CalculateDuplexMetrics.mean_family_size
      Float mean_ab_size = CalculateDuplexMetrics.mean_ab_size
      Float mean_ba_size = CalculateDuplexMetrics.mean_ba_size
      Float median_family_size = CalculateDuplexMetrics.median_family_size
      Float median_ab_size = CalculateDuplexMetrics.median_ab_size
      Float median_ba_size = CalculateDuplexMetrics.median_ba_size

      Float duplex_depth_above_500x = CollectStatisticsByCoverage.duplex_depth_above_500x
      Float duplex_depth_above_1000x = CollectStatisticsByCoverage.duplex_depth_above_1000x
      Float duplex_depth_above_1500x = CollectStatisticsByCoverage.duplex_depth_above_1500x
      Float duplex_depth_above_2000x = CollectStatisticsByCoverage.duplex_depth_above_2000x
   }
}


task GetBwaVersion {
   String bloodbiopsydocker
   String bwa_path
   Int? preemptible_attempts

   command {
      ${bwa_path} 2>&1 | \
      grep -e '^Version' | \
      sed 's/Version: //'
   }
   runtime {
      docker: bloodbiopsydocker
      memory: "1 GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 10])
   }
   output {
      String version = read_string(stdout())
   }
}

task DownsampleSam {
   String bloodbiopsydocker
   File? picard_override
   File bam_file
   File bam_index
   String base_name
   Float? downsample_probability
   Int? preemptible_attempts
   Int? memory
   Int disk_pad
   Int disk_size = ceil(size(bam_file, "GB") + size(bam_index, "GB")) + disk_pad
   Int mem = select_first([memory, 5]) 
   Int compute_mem = mem * 1000 - 500

   command {
      set -e

      export PICARD_LOCAL_JAR=${default="/usr/picard.jar" picard_override}

      if [ ${downsample_probability} == 1.0 ]; then
         cp ${bam_file} "${base_name}.bam"
         cp ${bam_index} "${base_name}.bai"
      else
         java -Xmx${compute_mem}m -jar $PICARD_LOCAL_JAR DownsampleSam \
            INPUT=${bam_file} \
            P=${downsample_probability} \
            OUTPUT=${base_name}.bam \
            CREATE_INDEX=true
      fi
   }
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 2])
   }
   output {
      File output_bam = "${base_name}.bam"
      File output_bam_index = "${base_name}.bai"
   }
}

task QuerySortSam {
   File input_bam
   File? picard_override
   String bloodbiopsydocker
   String base_name
   Int? preemptible_attempts
   Int? memory
   Int disk_pad
   Int disk_size = ceil(size(input_bam, "GB") * 5) + disk_pad
   Int mem = select_first([memory, 15])
   Int compute_mem = mem * 1000 - 500

  command {

   export PICARD_LOCAL_JAR=${default="/usr/picard.jar" picard_override}
      
   java -Xmx${compute_mem}m -jar $PICARD_LOCAL_JAR \
      SortSam \
      INPUT=${input_bam} \
      OUTPUT=${base_name}.sorted.bam \
      SORT_ORDER=queryname \
      MAX_RECORDS_IN_RAM=300000 \
      TMP_DIR=tmp 

  }
  output {
    File bam = "${base_name}.sorted.bam"
  }
  runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD, /cromwell_root/tmp 500 HDD"
      memory: mem + " GB"
  }
}

task AlignRawBamWithBwaMem {

   File input_bam
   String bwa_version
   String output_bam_basename
   File ref_fasta
   File ref_fasta_index
   File ref_dict

   # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
   # listing the reference contigs that are "alternative".
   File ref_alt

   File ref_amb
   File ref_ann
   File ref_bwt
   File ref_pac
   File ref_sa
   Int compression_level
   Int? preemptible_attempts
   Int? memory
   Int mem = select_first([memory, 15])
   Int disk_pad
   Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB"))
   Int disk_size = ceil(size(input_bam, "GB") * 5) + ref_size + disk_pad 

   command <<<
      set -o pipefail
      set -e

      # set the bash variable needed for the command-line
      bash_ref_fasta=${ref_fasta}
      # if ref_alt has data in it,
   
   
      java -Xms5000m -jar /usr/gitc/picard.jar \
          RevertSam \
          INPUT=${input_bam} \
          OUTPUT=reverted.bam
   
      java -Xms5000m -jar /usr/gitc/picard.jar \
          SamToFastq \
          INPUT=reverted.bam \
          FASTQ=/dev/stdout \
          INTERLEAVE=true \
          NON_PF=true | \
      /usr/gitc/bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta /dev/stdin - 2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) | \
      java -Dsamjdk.compression_level=${compression_level} -Xms3000m -jar /usr/gitc/picard.jar \
          MergeBamAlignment \
          VALIDATION_STRINGENCY=SILENT \
          EXPECTED_ORIENTATIONS=FR \
          ATTRIBUTES_TO_RETAIN=X0 \
          ATTRIBUTES_TO_REMOVE=NM \
          ATTRIBUTES_TO_REMOVE=MD \
          ALIGNED_BAM=/dev/stdin \
          UNMAPPED_BAM=reverted.bam \
          OUTPUT=${output_bam_basename}.unsorted.bam \
          REFERENCE_SEQUENCE=${ref_fasta} \
          PAIRED_RUN=true \
          SORT_ORDER="unsorted" \
          IS_BISULFITE_SEQUENCE=false \
          ALIGNED_READS_ONLY=false \
          CLIP_ADAPTERS=false \
          MAX_RECORDS_IN_RAM=2000000 \
          ADD_MATE_CIGAR=true \
          MAX_INSERTIONS_OR_DELETIONS=-1 \
          PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
          PROGRAM_RECORD_ID="bwamem" \
          PROGRAM_GROUP_VERSION="${bwa_version}" \
          PROGRAM_GROUP_COMMAND_LINE="bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta" \
          PROGRAM_GROUP_NAME="bwamem" \
          UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
          ALIGNER_PROPER_PAIR_FLAGS=true \
          UNMAP_CONTAMINANT_READS=true \
          ADD_PG_TAG_TO_READS=false 

      java -Xmx6000m -jar /usr/gitc/picard.jar \
         SortSam \
         INPUT=${output_bam_basename}.unsorted.bam \
         OUTPUT=${output_bam_basename}.bam \
         SORT_ORDER=coordinate \
         CREATE_INDEX=true


   >>>
   runtime {
      docker : "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135"
      preemptible: preemptible_attempts
      memory: mem + " GB"
      cpu: "16"
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
   }
   output {
      File output_bam = "${output_bam_basename}.bam"
      File output_bam_index = "${output_bam_basename}.bai"
      File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
   }
}

# This task takes a bam that has been deduplicated by
# start and stop positions and further deduplicates it
# by the UMI using the fgbio tool, GroupReadsByUmi.
# For more information about fgbio
# https://github.com/fulcrumgenomics/fgbio
task FGBioGroupReadsByUmi {
   String bloodbiopsydocker
   File? fgbio_override
   File bam_file
   File bam_index
   String base_name
   Int allowable_umi_distance
   Int? preemptible_attempts
   Int? memory
   Int disk_pad
   Int disk_size = ceil(size(bam_file, "GB") * 5) + ceil(size(bam_index, "GB")) + disk_pad 
   Int mem = select_first([memory, 5])
   Int compute_mem = mem * 1000 - 500

   command {
      set -e
      # It is necessary to set the tmp dir, otherwise
      # fgbio attempts to put the data in an invalid location.
      export FGBIO_LOCAL_JAR=${default="/usr/fgbio-1.0.0.jar" fgbio_override}

      java -Djava.io.tmpdir=/cromwell_root/ -Xmx${compute_mem}m -jar $FGBIO_LOCAL_JAR \
         GroupReadsByUmi \
         -e ${allowable_umi_distance} \
         -i ${bam_file} \
         -o ${base_name}.fgbio.groupByUmi.bam \
         --strategy=paired
   }
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 10])
   }
   output {
      File output_bam = "${base_name}.fgbio.groupByUmi.bam"
   }
}

# This task takes a bam that has been consensus called
# using fgbio CallDuplexConsensusReads.  
task FilterConsensusReads {
   String bloodbiopsydocker
   File? picard_override
   File? fgbio_override
   File reference
   String base_name
   File bam_file
   File bam_index
   String minimum_consensus_reads
   Int minimum_base_quality
   Float frac_Ns
   String? consensus_extra_filter_args
   Int? preemptible_attempts
   Int? memory
   Int disk_pad
   Int ref_size = ceil(size(reference, "GB"))
   Int disk_size = ceil(size(bam_file, "GB") * 5) + ceil(size(bam_index, "GB")) + ref_size + disk_pad 
   Int mem = select_first([memory, 5])
   Int compute_mem = mem * 1000 - 500

   command {
      set -e

      export FGBIO_LOCAL_JAR=${default="/usr/fgbio-1.0.0.jar" fgbio_override}
      export PICARD_LOCAL_JAR=${default="/usr/picard.jar" picard_override}

      # If the minimum consensus reads are set to 0, do nothing.
      if [ ${minimum_consensus_reads} == 0 ]; then
         cp ${bam_file} ${base_name}.fgbio.filterConsensus.bam
      else
         # At the time this was written FilterConsensusReads does not 
         # take --tmp-dir as an argument.
         # Because of this, we have to use -Djava.io.tmpdir to specify
         # a default tmp directory.  
         # Filters: 
         # - min-reads: the minimum number of reads supporting a consensus base
         # - min-base-quality: mask (make N) consensus bases with quality less than this threshold
         # - max-no-call-fraction: maximum fraction of no-calls in the read after filter
         java -Djava.io.tmpdir=/cromwell_root/tmp -Xmx${compute_mem}m -jar $FGBIO_LOCAL_JAR \
            FilterConsensusReads \
            -i ${bam_file} \
            -r ${reference} \
            -o ${base_name}.fgbio.filterConsensus.bam \
            --min-reads ${minimum_consensus_reads} \
            -n ${frac_Ns} \
            ${consensus_extra_filter_args} \
            --min-base-quality ${minimum_base_quality} 

      fi

      java -jar $PICARD_LOCAL_JAR BuildBamIndex \
         INPUT=${base_name}.fgbio.filterConsensus.bam \
         OUTPUT=${base_name}.fgbio.filterConsensus.bai

   }
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD, /cromwell_root/tmp 100 HDD"
      memory: mem + " GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 10])
   }
   output {
      File output_bam = "${base_name}.fgbio.filterConsensus.bam"
      File output_bam_index = "${base_name}.fgbio.filterConsensus.bai"
   }
}

#want to clip the end of the fragments 
#because they are lower quality 
#also want to filter for MQ < 60 
task ClipBam {
   String bloodbiopsydocker
   File? picard_override
   File? fgbio_override
   File bam_file
   File bam_index
   String base_name
   Int num_clip_bases_five_prime
   File reference
   File reference_index
   Int? num_clip_bases_three_prime
   Int? preemptible_attempts
   Int disk_pad 
   Int ref_size = ceil(size(reference, "GB") +  size(reference_index, "GB"))
   Int disk_size = ceil(size(bam_file, "GB") * 5) + ceil(size(bam_index, "GB")) + ref_size + disk_pad 
   Int? memory
   Int mem = select_first([memory, 10]) 
   Int compute_mem = mem * 1000 - 1000

   command {

      set -e 

      export FGBIO_LOCAL_JAR=${default="/usr/fgbio-1.0.0.jar" fgbio_override}
      export PICARD_LOCAL_JAR=${default="/usr/picard.jar" picard_override}

      java -Xmx${compute_mem}m -jar $FGBIO_LOCAL_JAR ClipBam \
         -i ${bam_file} \
         -o ${base_name}.clipped.bam \
         -c "Hard" \
         --ref ${reference} \
         --read-one-five-prime ${num_clip_bases_five_prime} \
         --read-two-five-prime ${num_clip_bases_five_prime} \
         ${"--read-one-three-prime " + num_clip_bases_three_prime} \
         ${"--read-two-three-prime " + num_clip_bases_three_prime} 

      samtools view -hb -q 60 ${base_name}.clipped.bam -o ${base_name}.filtered.bam 

      java -jar $PICARD_LOCAL_JAR BuildBamIndex \
         INPUT=${base_name}.filtered.bam  \
         OUTPUT=${base_name}.filtered.bai

   }
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 10])
   }
   output {
      File output_bam = "${base_name}.filtered.bam"
      File output_bam_index = "${base_name}.filtered.bai"
   }
}

# Collect HS or Targeted PCR Metrics
task CollectSelectionMetrics {
   String bloodbiopsydocker
   File? picard_override
   File reference
   File reference_index
   File reference_dict
   File bam_file
   File bam_index
   String base_name
   String suffix
   File bait_intervals
   File target_intervals
   Int? preemptible_attempts
   Int? memory
   Int disk_pad
   Int ref_size = ceil(size(reference, "GB") + size(reference_index, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(bam_file, "GB") * 2) + ceil(size(bam_index, "GB")) + ref_size + disk_pad 
   Int mem = select_first([memory, 5])
   Int compute_mem = mem * 1000 - 500

   command <<<
      set -e
      
      export PICARD_LOCAL_JAR=${default="/usr/picard.jar" picard_override}

      java -jar -Xmx${compute_mem}m $PICARD_LOCAL_JAR \
         CollectHsMetrics \
         I=${bam_file} \
         O=${base_name}${suffix}.selection_metrics \
         PER_TARGET_COVERAGE=${base_name}${suffix}.per_target_selection_metrics \
         THEORETICAL_SENSITIVITY_OUTPUT=${base_name}${suffix}.theoretical_sensitivity \
         $(for i in $(seq 0.001 .001 .02); do echo ALLELE_FRACTION=$i ;done) \
         BAIT_INTERVALS=${bait_intervals} \
         TARGET_INTERVALS=${target_intervals} \
         REFERENCE_SEQUENCE=${reference} \
         COVERAGE_CAP=1000000
      
   >>>
   runtime {
      docker: bloodbiopsydocker 
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 10])
   }
   output {
      File output_selection_metrics = "${base_name}${suffix}.selection_metrics"
      File output_per_target_selection_metrics = "${base_name}${suffix}.per_target_selection_metrics"
      File output_theoretical_sensitivity = "${base_name}${suffix}.theoretical_sensitivity"
   }
}


# Call consensus reads that have duplex UMI support.
task CallDuplexConsensusReads {
   String bloodbiopsydocker
   File bam_file
   File? fgbio_override
   String base_name
   String min_reads
   Int? preemptible_attempts
   Int? memory
   Int disk_pad
   Int disk_size = ceil(size(bam_file, "GB") * 5) + disk_pad 
   Int mem = select_first([memory, 10])
   Int compute_mem = mem * 1000 - 500

   command {
      set -e 

      export FGBIO_LOCAL_JAR=${default="/usr/fgbio-1.0.0.jar" fgbio_override}

      java -Xmx${compute_mem}m -Djava.io.tmpdir=/cromwell_root/tmp -jar $FGBIO_LOCAL_JAR \
         CallDuplexConsensusReads \
         -p ${base_name} \
         -i ${bam_file} \
         --min-reads ${min_reads}\
         -o ${base_name}.consensus.bam
   }
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD, /cromwell_root/tmp 100 HDD"
      memory: mem + " GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 10])
   }
   output {
      File output_bam = "${base_name}.consensus.bam"
   }
}


# Generates basic statistics on base quality score
# distribution.  Binning here is useful to reduce
# the numbers of categories used in recalibration.
# This task relies on R as well as Java which is why
# we use the genomes-in-the-cloud docker.
task QualityScoreDistribution {
   String bloodbiopsydocker
   File? picard_override
   File bam_file
   File bam_index
   String base_name
   Int? preemptible_attempts
   Int? memory
   Int disk_pad
   Int disk_size = ceil(size(bam_file, "GB") * 2.5) + ceil(size(bam_index, "GB")) + disk_pad 
   Int mem = select_first([memory, 5])
   Int compute_mem = mem * 1000 - 500

   command {

      set -e 

      export PICARD_LOCAL_JAR=${default="/usr/picard.jar" picard_override}

      java -Xmx${compute_mem}m -jar $PICARD_LOCAL_JAR QualityScoreDistribution \
         I=${bam_file} \
         O=${base_name}.qual_score_dist.txt \
         CHART=${base_name}.qual_score_dist.pdf
   }
   runtime {
      # We need this docker because the PDF is generated using an R-script
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 10])
   }
   output {
      File output_qual_score_dist = "${base_name}.qual_score_dist.txt"
      File output_qual_score_dist_pdf = "${base_name}.qual_score_dist.pdf"
   }
}

# Base quality score recalibration without the use of binned
# quality scores.  This is necessary after calling consensus
# reads because we want to know the precise quality scores.
task BQSRWithoutBinning {
   String bloodbiopsydocker
   File? gatk_override
   File bam_file
   File bam_index
   File reference
   File reference_index
   File reference_dict
   File dbsnp
   File dbsnp_index
   File known_indels
   File known_indels_index
   File variant_eval_gold_standard
   File variant_eval_gold_standard_index
   String base_name
   Int? preemptible_attempts
   Int? memory
   Int disk_pad
   Int ref_size = ceil(size(reference, "GB") + size(reference_index, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(bam_file, "GB") * 5) + ceil(size(bam_index, "GB")) + ref_size + disk_pad 
   Int mem = select_first([memory, 5])
   Int compute_mem = mem * 1000 - 500

   command {
      set -e

      export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

      gatk --java-options "-Xmx${compute_mem}m" BaseRecalibrator \
         -R ${reference} \
         -I ${bam_file} \
         --known-sites ${dbsnp} \
         --known-sites ${known_indels} \
         --known-sites ${variant_eval_gold_standard} \
         -O ${base_name}.recalibration_report.grp

      gatk --java-options "-Xmx${compute_mem}m" ApplyBQSR \
         -R ${reference} \
         -I ${bam_file} \
         -bqsr ${base_name}.recalibration_report.grp \
         -O ${base_name}.bqsr.bam
   }
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 10])
   }
   output {
      File output_bam = "${base_name}.bqsr.bam"
      File output_bam_index = "${base_name}.bqsr.bai"
      File output_recal_report_grp = "${base_name}.recalibration_report.grp"
   }
}

task AlignReadsAndMBA {

   String bloodbiopsydocker
   File? picard_override
   File reference
   File reference_bwa_index
   File reference_pac
   File reference_amb
   File reference_ann
   File reference_bwt
   File reference_sa
   File reference_dict
   File bam_file
   String base_name
   String bwa_path
   String bwa_version
   Int? preemptible_attempts
   Int? memory
   Int mem = select_first([memory, 32])
   Int disk_pad
   Int ref_size = ceil(size(reference, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(bam_file, "GB") * 15) + ref_size + disk_pad 

   command {
      set -e

      export PICARD_LOCAL_JAR=${default="/usr/picard.jar" picard_override}

      # Revert bam to fastq
      java -jar $PICARD_LOCAL_JAR SamToFastq \
         INPUT=${bam_file} \
         FASTQ=${base_name}.1.fastq \
         SECOND_END_FASTQ=${base_name}.2.fastq \
         UNPAIRED_FASTQ=${base_name}.unpaired.fastq

       # Align fastq with bwa
      ${bwa_path} mem -K 100000000 -t 16 \
         ${reference} \
         ${base_name}.1.fastq \
         ${base_name}.2.fastq > ${base_name}.aligned.sam

      java -jar -Xmx6144m $PICARD_LOCAL_JAR MergeBamAlignment \
         REFERENCE_SEQUENCE=${reference} \
         UNMAPPED_BAM=${bam_file} \
         ALIGNED_BAM=${base_name}.aligned.sam \
         OUTPUT=${base_name}.merged.bam \
         PROGRAM_RECORD_ID="bwamem" \
         PROGRAM_GROUP_VERSION="${bwa_version}" \
         PROGRAM_GROUP_COMMAND_LINE="bwa mem -K 100000000 -t 16" \
         PROGRAM_GROUP_NAME="bwamem" \
         CREATE_INDEX=true \
         TMP_DIR=.
   }
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + "GB"
      cpu: "16"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 10])
   }
   output {
      File output_bam = "${base_name}.merged.bam"
      File output_bam_index = "${base_name}.merged.bai"
   }
}



task CollectDepthOfCoverage {

   File interval_list 
   File reference
   File reference_index
   File reference_dict
   File bam_file
   File bam_index
   String base_name
   String? extra_arguments
   Int? preemptible_attempts
   Int? memory
   Int disk_pad
   Int ref_size = ceil(size(reference, "GB") + size(reference_index, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(bam_file, "GB") * 2) + ceil(size(bam_index, "GB")) + ref_size + disk_pad 
   Int mem = select_first([memory, 15])
   Int compute_mem = mem * 1000 - 1000

   command <<<
      set -e

      # Calculate tumor depth over the panel
      # only count fragments with same base.
      java -jar /usr/gitc/GATK36.jar -T DepthOfCoverage \
         -L ${interval_list} \
         -I ${bam_file} \
         -R ${reference} \
         -o ${base_name}.depth \
         --omitPerSampleStats \
         --printBaseCounts \
         ${extra_arguments}
   >>>

   runtime {
      docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135"
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + "GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 10])
   }

   output {
      File depth_of_coverage = "${base_name}.depth"
   }
}

task CollectStatisticsByCoverage {

   String bloodbiopsydocker
   String process_duplex_coverage_rscript 
   File raw_depth
   File start_stop_depth 
   File duplex_depth
   String base_name
   Int? preemptible_attempts
   Int? memory
   Int disk_pad
   Int disk_size = 200 
   Int mem = select_first([memory, 10]) 
   Int compute_mem = mem * 1000 - 500

   command <<<
      set -e

      Rscript -e 'source("${process_duplex_coverage_rscript}"); generateDepthFigures("${base_name}", "${raw_depth}", "${start_stop_depth}", "${duplex_depth}")'
   
      python <<CODE

      import pandas as pd

      df = pd.read_csv("${base_name}.depth.txt", delim_whitespace=True)

      def writeFile(value, filename):
         f = open(filename + ".txt", 'w')
         f.write(str(int(float(value))))
         f.close()

      for col in df.columns: 
         writeFile(df[col].loc[0], col)

      def writeDepthStatistic(depthValues, depthCutoff):
         f = open(str(depthCutoff) + ".txt", 'w')
         f.write(str(depthValues[depthValues.Total_Depth >= depthCutoff].count().Total_Depth / depthValues.count().Total_Depth))
         f.close()

      depthOfCoverageByLocus = pd.read_csv("${duplex_depth}", delimiter = '\t')

      writeDepthStatistic(depthOfCoverageByLocus, 500)
      writeDepthStatistic(depthOfCoverageByLocus, 1000)
      writeDepthStatistic(depthOfCoverageByLocus, 1500)
      writeDepthStatistic(depthOfCoverageByLocus, 2000)

      CODE

   >>>
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + "GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 10])
   }

   output {
      File depth_txt = "${base_name}.depth.txt"
      File raw_vs_duplex_depth = "${base_name}-RawVsDuplex.pdf"
      File start_stop_vs_raw_depth = "${base_name}-StartStopVsRawDepth.pdf"
      File start_stop_vs_duplex_depth = "${base_name}-StartStopVsDuplex.pdf"
      Int mean_raw_depth = read_int('meanRawDepth.txt')
      Int mean_startstop_depth = read_int('meanStartStopDepth.txt')
      Int mean_duplex_depth = read_int('meanDuplexDepth.txt')
      Float duplex_depth_above_500x = read_float('500.txt')
      Float duplex_depth_above_1000x = read_float('1000.txt')
      Float duplex_depth_above_1500x = read_float('1500.txt')
      Float duplex_depth_above_2000x = read_float('2000.txt')
   }

}

task CollectDuplexSeqMetrics {
   String bloodbiopsydocker
   File? fgbio_override
   File bam_file
   String base_name
   File target_intervals
   Int? preemptible_attempts
   Int? memory
   Int disk_pad
   Int disk_size = ceil(size(bam_file, "GB") * 2) + ceil(size(target_intervals, "GB"))  + disk_pad 
   Int mem = select_first([memory, 5])
   Int compute_mem = mem * 1000 - 500
   
   command {

      set -e 

      export FGBIO_LOCAL_JAR=${default="/usr/fgbio-1.0.0.jar" fgbio_override}

      java -Xmx${compute_mem}m -jar $FGBIO_LOCAL_JAR CollectDuplexSeqMetrics \
         -i ${bam_file} \
         -l ${target_intervals} \
         --output ${base_name}
   }
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 10])
   }
   output {
      File duplex_family_sizes = "${base_name}.duplex_family_sizes.txt"
      File duplex_qc = "${base_name}.duplex_qc.pdf"
      File duplex_yield_metrics = "${base_name}.duplex_yield_metrics.txt"
      File family_sizes = "${base_name}.family_sizes.txt"
      File umi_counts = "${base_name}.umi_counts.txt"
   }
}

task CalculateDuplexMetrics {

   String bloodbiopsydocker
   String calcFS_script
   File duplex_family_sizes
   Int disk_size = 200
   Int? preemptible_attempts
   Int? memory
   Int mem = select_first([memory, 5])

   command <<<

      Rscript -e "source('${calcFS_script}');getDuplexFamilySizeMetrics('${duplex_family_sizes}')"

   >>>
   runtime {
      docker: bloodbiopsydocker
      memory: mem + " GB"
      disks: "local-disk " + disk_size + " HDD"
      maxRetries: 3
   }
   output {
      Float mean_family_size = read_float("mean_total_family_size.txt")
      Float mean_ab_size = read_float("mean_ab_family_size.txt")
      Float mean_ba_size = read_float("mean_ba_family_size.txt")
      Float median_family_size = read_float("median_total_family_size.txt")
      Float median_ab_size = read_float("median_ab_family_size.txt")
      Float median_ba_size = read_float("median_ba_family_size.txt")

   }
}
