version 1.0

workflow Interval2Bed {
    input {
        Array[File] interval_files
        Array[String]? interval_labels
    }

    if (defined(interval_labels)) {
        Array[Pair[String, File]] paired_intervals = zip(select_first([interval_labels]), interval_files)

        scatter (p in paired_intervals) {
            # Check that file extension is .interval_list, i.e. removing it shortens the basename
            if (basename(p.right, ".interval_list") != basename(p.right)) {
                call Convert as ConvertPairs {
                    input:
                        interval_list=p.right,
                        interval_label=p.left
                }
            }

            if (basename(p.right, ".bed") != basename(p.right)) {
                File bed_file_pair = p.right
                File bed_label_pair = p.left
            }
        }

        Array[File] new_interval_files_pair = flatten([select_all(bed_file_pair), select_all(ConvertPairs.bed_file)])
        Array[File] new_interval_labels_pair = flatten([select_all(bed_label_pair), select_all(ConvertPairs.bed_label)])
    }

    if (!defined(interval_labels)) {
        scatter (f in interval_files) {
            if (basename(f, ".interval_list") != basename(f)) {
                call Convert as ConvertUnlabeled {
                    input:
                        interval_list=f,
                        interval_label=basename(f, ".interval_list")
                }
            }

            if (basename(f, ".bed") != basename(f)) {
                File bed_file = f
                String bed_label = basename(f, ".bed")
            }
        }

        Array[File] new_interval_files = flatten([select_all(bed_file), select_all(ConvertUnlabeled.bed_file)])
        Array[File] new_interval_labels = flatten([select_all(bed_label), select_all(ConvertUnlabeled.bed_label)])
    }


    output {
        Array[File] bed_files = select_first([new_interval_files_pair, new_interval_files])
        Array[String] bed_labels = select_first([new_interval_labels_pair, new_interval_labels])
    }
}

task Convert {
    input {
        File interval_list
        String? interval_label

        String gatk_tag = "4.2.3.0"
        Int? preemptible
        Int disk_size = ceil(4 * size(interval_list, "GB")) + 5
        Int cpu = 4
        Int memory = 16
    }

    String name = basename(interval_list, ".interval_list")

    command <<<
        set -xe

        gatk IntervalListToBed -I ~{interval_list} -O ~{name}.bed
    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:" + gatk_tag
        preemptible: select_first([preemptible, 0])
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: "16"
        cpu: cpu
        memory: memory + " GB"
    }

    output {
        File bed_file = "~{name}.bed"
        String? bed_label = select_first([interval_label, name])
    }
}