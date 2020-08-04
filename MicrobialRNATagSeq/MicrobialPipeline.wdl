version 1.0

import "tasks.wdl" as tasks
import "structs.wdl"

workflow MicrobialPipeline {
    input {
        Array[File] fastq_r1s
        Array[File] fastq_r2s
        Array[String]+ sample_names
        Array[String]+ barcodes
        Array[String]+ libraries
        String project_id
        Array[File]+ gffs
        Array[File]+ fastas
        Array[String]+ microbeNames
        Array[String] types = ["CDS", "ncRNA", "rRNA", "tRNA"]
    }

    call tasks.GetBwaVersion

    call tasks.BuildCombinedName {
        input:
            names = microbeNames
    }

    call tasks.BuildReferenceResources {
        input:
            fastas = fastas,
            gffs = gffs,
            combined_name = BuildCombinedName.combined_name

    }

    call tasks.BuildSampleSheet {
        input:
            sample_names = sample_names,
            barcodes = barcodes,
            libraries = libraries
    }

    ReferenceFasta ref = object {
                        ref_dict : BuildReferenceResources.combined_dict,
                        ref_fasta : BuildReferenceResources.combined_fasta,
                        ref_fasta_index : BuildReferenceResources.combined_index,
                        ref_sa : BuildReferenceResources.combined_sa,
                        ref_amb : BuildReferenceResources.combined_amb,
                        ref_bwt : BuildReferenceResources.combined_bwt,
                        ref_ann : BuildReferenceResources.combined_ann,
                        ref_pac : BuildReferenceResources.combined_pac,
                                }


    scatter (i in range(length(fastq_r1s))) {
        FastQSet fastq_set = object {
            fastq_r1 : fastq_r1s[i],
            fastq_r2 : fastq_r2s[i]
        }
        call tasks.SplitBarcodes {
            input:
                fastq_set = fastq_set,
                sample_sheet = BuildSampleSheet.sample_sheet
        }

        scatter (unmapped_bam in SplitBarcodes.sample_unmapped_bams) {
            call tasks.AlignPairedReadsBWABacktrack {
                input:
                    unmapped_bam = unmapped_bam,
                    ref = ref,
                    bwa_version = GetBwaVersion.bwa_version
            }
        }
    }

    Array[Array[File]] bams_to_merge_arrays = transpose(AlignPairedReadsBWABacktrack.aligned_bam)
    Array[Array[File]] bais_to_merge_arrays = transpose(AlignPairedReadsBWABacktrack.aligned_bam_index)

    scatter (i in range(length(bams_to_merge_arrays))) {
        Array[File] bams_to_merge = bams_to_merge_arrays[i]
        Array[File] bais_to_merge = bais_to_merge_arrays[i]

        call tasks.GetSampleName {
            input:
                bam = bams_to_merge[0]
        }
        call tasks.MergeBamsByBarcodes {
            input:
                aligned_bams = bams_to_merge,
                aligned_bam_indecies = bais_to_merge,
                sample_name = GetSampleName.sample_name
        }

        call tasks.CountFeaturesGATK {
            input:
                aligned_bam = MergeBamsByBarcodes.merged_bam,
                aligned_bam_index = MergeBamsByBarcodes.merged_bam_index,
                feature_file = BuildReferenceResources.combined_gff,
                gatk_jar = "gs://broad-dsde-methods-ckachulis/jars/gatk_collect_fragment_counts_microbial_v5.jar",
                sample_name = GetSampleName.sample_name,
                sequence_dict = ref.ref_dict,
                types=types
        }

        call tasks.CountsToFPKMGATK {
                    input:
                        counts = CountFeaturesGATK.counts,
                        alignmentMetrics = CollectAlignmentSummaryMetrics.alignment_summary_metrics
        }

        call tasks.BamToTDF {
            input:
                bam = MergeBamsByBarcodes.merged_bam,
                bam_index = MergeBamsByBarcodes.merged_bam_index
        }

        call tasks.CollectAlignmentSummaryMetrics {
            input:
                bam = MergeBamsByBarcodes.merged_bam,
                bam_index = MergeBamsByBarcodes.merged_bam_index,
                ref_fasta = ref.ref_fasta
        }

        call tasks.MeanQualityByCycle {
            input:
                bam = MergeBamsByBarcodes.merged_bam,
                bam_index = MergeBamsByBarcodes.merged_bam_index
        }

        call tasks.CollectInsertSizeMetrics {
            input:
                bam = MergeBamsByBarcodes.merged_bam,
                bam_index = MergeBamsByBarcodes.merged_bam_index
        }

        call tasks.CollectRnaSeqMetrics {
            input:
                bam = MergeBamsByBarcodes.merged_bam,
                bam_index = MergeBamsByBarcodes.merged_bam_index,
                gff = BuildReferenceResources.combined_gff,
                picardJar = "gs://broad-dsde-methods-ckachulis/jars/picard_CollectRNAMetrics.jar"
        }
    }

    call tasks.CombineGeneCountsGATK as CombineCounts {
        input:
            count_files = CountFeaturesGATK.counts,
            outName = "All.counts"
    }

    call tasks.CombineGeneCountsGATK as CombineFPKM {
        input:
            count_files = CountsToFPKMGATK.fpkm,
            outName = "All.fpkm"
    }

    call tasks.SampleCorrelationsGATK {
        input:
            fpkm = CombineFPKM.combined
    }

    call tasks.CombineMetrics as CombineAlignmentSummaryMetrics {
        input:
            metrics_files = CollectAlignmentSummaryMetrics.alignment_summary_metrics
    }

    call tasks.CombineMetrics as CombineInsertSizeMetrics {
        input:
            metrics_files = CollectInsertSizeMetrics.insert_size_metrics
    }

    call tasks.CombineMetrics as CombineRnaSeqMetrics {
        input:
            metrics_files = CollectRnaSeqMetrics.rna_seq_metrics
    }

    output {
        File rna_seq_metrics = CombineRnaSeqMetrics.combined_metrics
        File alignment_summary_metrics = CombineAlignmentSummaryMetrics.combined_metrics
        File insert_size_metrics = CombineInsertSizeMetrics.combined_metrics
        File combined_counts = CombineCounts.combined
        File combined_fpkm = CombineFPKM.combined
        File correlations = SampleCorrelationsGATK.corr_file
    }

#    call tasks.Picard_metrics_parse {
#        input:
#            metrics = CollectAlignmentSummaryMetrics.alignment_summary_metrics,
#            project = project_id
#    }

#    call tasks.RPG_metrics {
#        input:
#            counts = CountFeaturesJL.counts,
#            mets = CountFeaturesJL.mets,
#            project_id = project_id
#    }
}






