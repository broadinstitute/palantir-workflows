version 1.0

task CalculateCoverage {
   input {
       File bam_file
       File bam_index
       File interval_list
       String base_name
       Int min_base_quality
       Int min_mapping_quality
       File ref_fasta
       File ref_dict
       File ref_index
       String? java_mem
   }

    Int disk_size = 10 + ceil(size(bam_file, "GiB") + 10*size(ref_fasta, "GB") + size(ref_index,"GB"))
   command <<<
      java ~{"-Xmx"+java_mem} -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage \
         -L ~{interval_list} \
         -I ~{bam_file} \
         -mmq ~{min_mapping_quality} \
         -mbq ~{min_base_quality} \
         --countType COUNT_FRAGMENTS \
         -o "~{base_name}.coverage" \
         -R ~{ref_fasta}
   >>>
   output {
      File coverage = "~{base_name}.coverage"
   }
   runtime {
      docker: "broadinstitute/gatk3:3.8-1"
      disks: "local-disk " + disk_size + " HDD"
      memory: "16 GB"
   }
}

task GetSampleName {
    input {
        File bam
    }

    Int disk_size = 10 + ceil(size(bam, "GiB"))
    command <<<
        set -e

        samtools view -H ~{bam} | \
            awk '$1=="@RG" {
                    for(i=2; i<=NF; i++) {
                        split($i,split_str,":")
                        if(split_str[1]=="SM") {
                            print split_str[2]
                        }
                    }
                }' | sort | uniq
    >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
    disks: "local-disk "+disk_size+" HDD"
    memory: "1 GB"
  }

  output {
    String sample_name = read_string(stdout())
  }
}

task CollectData {
    input {
       Array[File] coverage_files
       File sexes_file
       String base_name
       Int low_coverage_threshold
       Float sample_threshold
   }

	Int disk_size = 10 + ceil(1.1 * size(coverage_files, "GiB"))
   command <<<
      set -e
      set -v

      # Paste all the coverage files together into a single file base_name.coverage.txt
      echo 'paste <(cat ~{coverage_files[0]} | cut -f 1) \' > run.sh

      echo ~{sep=" " coverage_files} > coverage_files.txt
      sed -i 's/ /\n/g' coverage_files.txt

      while read line
      do
         echo '<(cat' $line '| cut -f 4) \' >> run.sh
      done < coverage_files.txt
      chmod 755 run.sh
      ./run.sh > ~{base_name}.coverage.txt

      echo "Starting Data Munging in R"

      # In R do various data munging operations
      Rscript -<<"EOF"
      library(data.table)
      sites <- fread("~{base_name}.coverage.txt")
      sexes <- fread("~{sexes_file}", col.names = c("sample", "sex", "normalizedX", "normalizedY"))
      maleColumns <- paste0("Depth_for_", sexes[sex == "Male", sample])
      femaleColumns <- paste0("Depth_for_", sexes[sex == "Female", sample])

      sites[, contig := gsub(":.*", "", Locus)]
      sites[, pos := as.numeric(gsub(".*:", "", Locus))]
      sites[, countBelowThreshold := Reduce(`+`, lapply(.SD, "<", ~{low_coverage_threshold})), .SDcols = c(maleColumns, femaleColumns)]
      sites[grepl("^(chr)?X", contig), maleCountBelowThreshold := Reduce(`+`, lapply(.SD, "<", ~{low_coverage_threshold} / 2)), .SDcols = maleColumns]
      sites[grepl("^(chr)?X", contig), femaleCountBelowThreshold := Reduce(`+`, lapply(.SD, "<", ~{low_coverage_threshold})), .SDcols = femaleColumns]
      sites[grepl("^(chr)?X", contig), countBelowThreshold := femaleCountBelowThreshold + maleCountBelowThreshold]
      sites[grepl("^(chr)?Y", contig), countBelowThreshold := Reduce(`+`, lapply(.SD, "<", ~{low_coverage_threshold} / 2)), .SDcols = maleColumns]

      sampleThreshold <- sexes[, .N] * ~{sample_threshold}
      maleSampleThreshold <- sexes[sex == "Male", .N] * ~{sample_threshold}

      sites[!grepl("^(chr)?Y", contig), poorlyCovered := countBelowThreshold >= sampleThreshold]
      sites[grepl("^(chr)?Y", contig), poorlyCovered := countBelowThreshold >= maleSampleThreshold]
      poorlyCoveredSites <- sites[poorlyCovered == TRUE, .(contig, pos)]
      poorlyCoveredSites[, interval := cumsum(c(1, (pos - shift(pos) > 1 | contig != shift(contig))[-1]))]
      poorlyCoveredIntervals <- unique(poorlyCoveredSites[, .(contig = contig, start = min(pos) - 1, end = max(pos)), by = interval])
      poorlyCoveredIntervals[, orientation := "+"]
      poorlyCoveredIntervals[, name := paste0(contig, "_", start, "_", end)]
      poorlyCoveredIntervals[, interval := NULL]
      poorlyCoveredIntervals
      fwrite(poorlyCoveredIntervals, "~{base_name}.low_coverage.bed", col.names = FALSE, sep = "\t")
      EOF

   >>>
   output {
      File coverage_output = "~{base_name}.coverage.txt"
      File bed_file_output = "~{base_name}.low_coverage.bed"
   }
   runtime {
      docker: "rocker/tidyverse:3.5.2"
      disks: "local-disk " + disk_size + " HDD"
      memory: "100 GB"
   }
}

task DetermineXYCoverage {
    input {
        File bam_file
        File bam_index
        String base_name
        Int disk_size
    }

    command <<<
          samtools idxstats ~{bam_file} > ~{base_name}.idxstats.txt

          Rscript -<<"EOF"
          idxstats <- read.table("~{base_name}.idxstats.txt", col.names = c("contig", "length", "mapped", "unmapped"))

          # Determine normalization constant.  This is the number of bases/chromosome/read.
          normalization <- sum(as.numeric(idxstats[grepl("^(chr)?[0-9]{1,2}$", idxstats$contig), ]$length)) /
              sum(as.numeric(idxstats[grepl("^(chr)?[0-9]{1,2}$", idxstats$contig), ]$mapped))
          # Estimate the number of X and Y chromosomes in the sample by
          # normalizing the coverage on X and Y.
          normalizedX <- 2 / (idxstats[grepl("^(chr)?X$", idxstats$contig), ]$length /
                                  idxstats[grepl("^(chr)?X$", idxstats$contig), ]$mapped / normalization)
          normalizedY <- 2 / (idxstats[grepl("^(chr)?Y$", idxstats$contig), ]$length /
                                  idxstats[grepl("^(chr)?Y$", idxstats$contig), ]$mapped / normalization)
          sampleAndXYChroms <- data.frame("~{base_name}", round(normalizedX, 3), round(normalizedY, 3))

          write.table(sampleAndXYChroms, "~{base_name}.xy.txt", col.name = F, row.name = F, quote = F)
          EOF
    >>>

    output {
          File xy_out = "~{base_name}.xy.txt"
    }

    runtime {
          docker: "fleharty/genomes@sha256:d162bf0de0b470ec3bdd2ffc307186a27e47e4c001b6f21e901515d2451dd694"
          disks: "local-disk " + disk_size + " HDD"
          memory: "16 GB"
       }
}

task DetermineSexesByClustering {
    input {
        Array[File] inputFiles
    }

    command <<<
    Rscript -<<"EOF"
    library(data.table)
    d_list <- lapply(list("~{sep="\",\"" inputFiles}"), fread)
    d <- rbindlist(d_list)
    d_mat <- as.matrix(d,rownames = "V1")
    initial_centers <- matrix(c(1,2,1,0), nrow=2)
    clusters <- kmeans(d_mat,initial_centers)
    d[,cluster:=list(clusters$cluster)]
    d_centers <- data.table(clusters$centers)
    d_centers[,d_male:=sqrt((V2-1)^2 + (V3-1)^2)]
    d_centers[,d_female:=sqrt((V2-2)^2 + (V3-0)^2)]
    if (d_centers[1,d_male]<d_centers[2,d_male] && d_centers[1,d_female]>d_centers[2,d_female]) {
        d[cluster==2,sex:="Female"]
        d[cluster==1,sex:="Male"]
    } else if (d_centers[1,d_male]>d_centers[2,d_male] && d_centers[1,d_female]<d_centers[2,d_female]) {
        d[cluster==2,sex:="Male"]
        d[cluster==1,sex:="Female"]
    } else {
        d[,sex:="Unknown"]
    }
    d[,cluster:=NULL]
    setcolorder(d,c("V1","sex","V2","V3"))
    write.table(d, "all.sexes.txt", col.name = F, row.name = F, quote = F)
    EOF

    >>>

    runtime {
        docker: "rocker/tidyverse:3.5.2"
    }

    output {
        File output_sexes = "all.sexes.txt"
    }
}

task BedToIntervalList {
    input {
        File bed
        File sequence_dict
    }

    String output_name = basename(bed, ".bed") + ".interval_list"
    Int disk_size = 10 + ceil(2.5 * size(bed, "GB"))

    command <<<
        java -jar /usr/gitc/picard.jar BedToIntervalList I=~{bed} SD=~{sequence_dict} O=~{output_name}
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
        disks: "local-disk "+disk_size+" HDD"
    }

    output {
        File interval_list = "~{output_name}"
    }
}

workflow GetUndercoveredRegions {
    input {
       Array[File] bam_files
       Array[File] bam_indicies
       String base_name
       File interval_list
       Int min_base_quality
       Int min_mapping_quality
       Int low_coverage_threshold
       Float sample_fraction_threshold
       Int disk_size
       File ref_dict
       File ref_fasta
       File ref_index
       String? java_mem
   }

   Array[Pair[File, File]] bam_index_pairs = zip(bam_files, bam_indicies)


   scatter (paired_data in bam_index_pairs) {

      String basename = basename(paired_data.left, ".cram")

      call GetSampleName {
        input:
            bam = paired_data.left
      }

      call DetermineXYCoverage {
         input:
            bam_file = paired_data.left,
            bam_index = paired_data.right,
            base_name = GetSampleName.sample_name,
            disk_size = 200
      }

      call CalculateCoverage {
         input:
            bam_file = paired_data.left,
            bam_index = paired_data.right,
            base_name = basename,
            interval_list = interval_list,
            min_base_quality = min_base_quality,
            min_mapping_quality = min_mapping_quality,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            java_mem = java_mem
      }
   }

   call DetermineSexesByClustering {
        input:
            inputFiles = DetermineXYCoverage.xy_out
   }

   call CollectData {
      input:
        coverage_files = CalculateCoverage.coverage,
        sexes_file = DetermineSexesByClustering.output_sexes,
        base_name = base_name,
        low_coverage_threshold = low_coverage_threshold,
        sample_threshold = sample_fraction_threshold
   }

   call BedToIntervalList {
      input:
        bed = CollectData.bed_file_output,
        sequence_dict = ref_dict
   }

   output {
      File coverage_output = CollectData.coverage_output
      File bed_file_output = CollectData.bed_file_output
      File interval_list_output = BedToIntervalList.interval_list
      File all_sexes_output = DetermineSexesByClustering.output_sexes
   }
}

