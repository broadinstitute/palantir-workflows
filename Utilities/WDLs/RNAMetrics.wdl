version 1.0

import "RNAMetricsTasks.wdl" as RNAMetricsTasks

workflow RNAMetrics {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceAnnotation
        File referenceGenome
        File referenceGenomeIndex
        File refFlat
        File ribosomalIntervals
    }

    call RNAMetricsTasks.RNASeQC2 as RNASeQC2 {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceAnnotation = referenceAnnotation
    }

    call RNAMetricsTasks.CollectRNASeqMetrics as CollectRNASeqMetrics {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            refFlat = refFlat,
            ribosomalIntervals = ribosomalIntervals
    }

    call RNAMetricsTasks.MultiQC as MultiQC {
        input:
            collectRNASeqMetricsOutput = CollectRNASeqMetrics.rnaMetrics,
            rnaSeQCOutput = RNASeQC2.metrics
    }

    output {
        File exonCV = RNASeQC2.exonCV
        File exonReads = RNASeQC2.exonReads
        File geneFragments = RNASeQC2.geneFragments
        File geneReads = RNASeQC2.geneReads
        File geneTPM = RNASeQC2.geneTPM
        File metrics = RNASeQC2.metrics
        File rnaseqcMonitoringLog = RNASeQC2.monitoringLog
        File rnaMetrics = CollectRNASeqMetrics.rnaMetrics
        File collectRNASeqMetricsMonitoringLog = CollectRNASeqMetrics.monitoringLog
        File multiQCReport = MultiQC.multiQCReport
        File multiQCData = MultiQC.multiQCData
    }
}