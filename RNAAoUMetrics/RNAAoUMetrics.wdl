version 1.0

import "RNAAoUMetricsTasks.wdl" as RNAAoUMetricsTasks

workflow RNAAoUMetrics {
    input {
        File alignment
        File alignmentIndex
        File reference
        File referenceIndex
        File refFlat
        File ribosomalIntervals
    }

    call RNAAoUMetricsTasks.CollectRNASeqMetrics {
        input:
            alignment = alignment,
            alignmentIndex = alignmentIndex,
            reference = reference,
            referenceIndex = referenceIndex,
            refFlat = refFlat,
            ribosomalIntervals = ribosomalIntervals
    }

    output {
        File rnaSeqMetrics = CollectRNASeqMetrics.rnaSeqMetrics
    }
}