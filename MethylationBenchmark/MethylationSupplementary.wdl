version 1.0

import "MethylationSupplementaryTasks.wdl" as MethylationSupplementaryTasks

workflow MethylationSupplementary {
    input {
        String sampleId
        File bedGraph
        File chromSizes
    }

    call MethylationSupplementaryTasks.CreateBigWigFromBedGraph {
        input:
            sampleId = sampleId,
            bedGraph = bedGraph,
            chromSizes = chromSizes
    }

    call MethylationSupplementaryTasks.MultiBigWigSummary {

    }

    output {
        File bigWig = CreateBigWigFromBedGraph.bigWig
        File multiBigWigSummary = MultiBigWigSummary.multiBigWigSummary
    }
}