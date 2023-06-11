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

    call RNAAoUMetricsTasks.CollectInsertSizeMetrics {
        input:
            alignment = alignment,
            alignmentIndex = alignmentIndex
    }

    call RNAAoUMetricsTasks.CollectAlignmentSummaryMetrics {
        input:
            alignment = alignment,
            alignmentIndex = alignmentIndex,
            reference = reference,
            referenceIndex = referenceIndex
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
        File insertSizeMetrics = CollectInsertSizeMetrics.insertSizeMetrics
        File insertSizeHistogram = CollectInsertSizeMetrics.insertSizeHistogram
        File alignmentSummaryMetrics = CollectAlignmentSummaryMetrics.alignmentSummaryMetrics
        File rnaSeqMetrics = CollectRNASeqMetrics.rnaSeqMetrics
    }
}