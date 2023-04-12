version 1.0

import "LongReadRNABenchmarkTasks.wdl" as LongReadRNABenchmarkTasks

workflow LongReadRNABenchmark {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        File expressedGTF
        File expressedKeptGTF
        File excludedGTF
        String datasetName
        String dataType
    }

    call LongReadRNABenchmarkTasks.IsoQuant as IsoQuant {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            dataType = dataType,
            cpu = 16,
            numThreads = 32,
            memoryGB = 256,
            diskSizeGB = 500,
            docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant:latest",
            monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

#    call LongReadRNABenchmarkTasks.IsoQuant as IsoQuantReferenceFree {
#        input:
#            inputBAM = inputBAM,
#            inputBAMIndex = inputBAMIndex,
#            referenceGenome = referenceGenome,
#            referenceGenomeIndex = referenceGenomeIndex,
#            datasetName = datasetName,
#            dataType = dataType,
#            cpu = 16,
#            numThreads = 32,
#            memoryGB = 256,
#            diskSizeGB = 500,
#            docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant:latest",
#            monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
#    }

    call LongReadRNABenchmarkTasks.StringTie as StringTie {
        input:
            inputBAM = inputBAM,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            cpu = 16,
            numThreads = 32,
            memoryGB = 64,
            diskSizeGB = 500,
            docker = "us.gcr.io/broad-dsde-methods/kockan/stringtie:latest",
            monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

#    call LongReadRNABenchmarkTasks.StringTie as StringTieReferenceFree {
#        input:
#            inputBAM = inputBAM,
#            datasetName = datasetName,
#            cpu = 16,
#            numThreads = 32,
#            memoryGB = 64,
#            diskSizeGB = 500,
#            docker = "us.gcr.io/broad-dsde-methods/kockan/stringtie:latest",
#            monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
#    }

    call LongReadRNABenchmarkTasks.Bambu as Bambu {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            dataType = dataType,
            cpu = 16,
            numThreads = 32,
            memoryGB = 64,
            diskSizeGB = 500,
            docker = "us.gcr.io/broad-dsde-methods/kockan/bambu:latest",
            monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    call LongReadRNABenchmarkTasks.Flair as Flair {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            cpu = 16,
            numThreads = 32,
            memoryGB = 64,
            diskSizeGB = 500,
            docker = "brookslab/flair:latest",
            monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    call LongReadRNABenchmarkTasks.Talon as Talon {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            dataType = dataType,
            cpu = 16,
            numThreads = 32,
            memoryGB = 256,
            diskSizeGB = 500,
            docker = "us.gcr.io/broad-dsde-methods/kockan/talon:latest",
            monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    call LongReadRNABenchmarkTasks.IsoSeq as IsoSeq {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            datasetName = datasetName,
            cpu = 16,
            numThreads = 32,
            memoryGB = 256,
            diskSizeGB = 500,
            docker = "us.gcr.io/broad-dsde-methods/kockan/isoseq3:latest",
            monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    call LongReadRNABenchmarkTasks.Tama as Tama {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            datasetName = datasetName,
            cpu = 16,
            memoryGB = 256,
            diskSizeGB = 500,
            docker = "us.gcr.io/broad-dsde-methods/kockan/tama:latest",
            monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    call LongReadRNABenchmarkTasks.Flames as Flames {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            cpu = 16,
            memoryGB = 256,
            diskSizeGB = 500,
            docker = "us.gcr.io/broad-dsde-methods/kockan/flames:latest",
            monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    call LongReadRNABenchmarkTasks.Cupcake as Cupcake {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            datasetName = datasetName,
            cpu = 16,
            memoryGB = 256,
            diskSizeGB = 500,
            docker = "us.gcr.io/broad-dsde-methods/kockan/cdna-cupcake:latest",
            monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareIsoQuant {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF,
            inputGTF = IsoQuant.isoQuantGTF,
            toolName = "isoquant",
            datasetName = datasetName,
            cpu = 8,
            memoryGB = 64,
            diskSizeGB = 300,
            docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareStringTie {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF,
            inputGTF = StringTie.stringTieGTF,
            toolName = "stringtie",
            datasetName = datasetName,
            cpu = 8,
            memoryGB = 64,
            diskSizeGB = 300,
            docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareBambu {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF,
            counts = Bambu.bambuGTFCounts,
            inputGTF = Bambu.bambuGTF,
            toolName = "bambu",
            datasetName = datasetName,
            cpu = 8,
            memoryGB = 64,
            diskSizeGB = 300,
            docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareFlair {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF,
            inputGTF = Flair.flairGTF,
            toolName = "flair",
            datasetName = datasetName,
            cpu = 8,
            memoryGB = 64,
            diskSizeGB = 300,
            docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareTalon {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF,
            inputGTF = Talon.talonGTF,
            toolName = "talon",
            datasetName = datasetName,
            cpu = 8,
            memoryGB = 64,
            diskSizeGB = 300,
            docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
    }

#    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareIsoSeq {
#        input:
#            reducedAnnotationDB = IsoQuant.isoQuantDB,
#            expressedGTF = expressedGTF,
#            expressedKeptGTF = expressedKeptGTF,
#            excludedGTF = excludedGTF,
#            inputGTF = IsoSeq.isoSeqGFF,
#            toolName = "isoseq",
#            datasetName = datasetName,
#            cpu = 8,
#            memoryGB = 64,
#            diskSizeGB = 300,
#            docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
#    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareFlames {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF,
            inputGTF = Flames.flamesGFF,
            toolName = "flames",
            datasetName = datasetName,
            cpu = 8,
            memoryGB = 64,
            diskSizeGB = 300,
            docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
    }

#    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareCupcake {
#        input:
#            reducedAnnotationDB = IsoQuant.isoQuantDB,
#            expressedGTF = expressedGTF,
#            expressedKeptGTF = expressedKeptGTF,
#            excludedGTF = excludedGTF,
#            inputGTF = Cupcake.cupcakeGFF,
#            toolName = "cupcake",
#            datasetName = datasetName,
#            cpu = 8,
#            memoryGB = 64,
#            diskSizeGB = 300,
#            docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
#    }

#    call LongReadRNABenchmarkTasks.DenovoAnnotationGFFCompare as DenovoAnnotationGFFCompare {
#        input:
#            isoQuantGTF = IsoQuant.isoQuantGTF,
#            stringTieGTF = StringTie.stringTieGTF,
#            bambuGTF = Bambu.bambuGTF,
#            bambuGTFCounts = Bambu.bambuGTFCounts,
#            flairGTF = Flair.flairGTF,
#            talonGTF = Talon.talonGTF,
#            datasetName = datasetName
#    }

#    call LongReadRNABenchmarkTasks.ReferenceFreeGFFCompare as ReferenceFreeGFFCompare {
#        input:
#            isoQuantDenovoGTF = IsoQuantReferenceFree.isoQuantDenovoGTF,
#            stringTieDenovoGTF = StringTieReferenceFree.stringTieDenovoGTF,
#            expressedGTF = expressedGTF,
#            datasetName = datasetName
#    }

    call LongReadRNABenchmarkTasks.ReducedAnalysisSummarize as ReducedAnalysisSummarize {
        input:
            reducedGffCompareOutIsoQuant = ReducedAnnotationGFFCompareIsoQuant.gffCompareOutput,
            reducedGffCompareOutStringTie = ReducedAnnotationGFFCompareStringTie.gffCompareOutput,
            reducedGffCompareOutBambu = ReducedAnnotationGFFCompareBambu.gffCompareOutput,
            reducedGffCompareOutFlair = ReducedAnnotationGFFCompareFlair.gffCompareOutput,
            reducedGffCompareOutTalon = ReducedAnnotationGFFCompareTalon.gffCompareOutput,
#            reducedGffCompareOutIsoSeq = ReducedAnnotationGFFCompareIsoSeq.gffCompareOutput,
            reducedGffCompareOutFlames = ReducedAnnotationGFFCompareFlames.gffCompareOutput,
#            reducedGffCompareOutCupcake = ReducedAnnotationGFFCompareCupcake.gffCompareOutput,
            datasetName = datasetName,
            cpu = 1,
            memoryGB = 32,
            diskSizeGB = 100,
            docker = "us.gcr.io/broad-dsde-methods/kockan/kockan-reduced-analysis-summarize:latest"
    }

    output {
        File isoQuantDB = IsoQuant.isoQuantDB
        File isoQuantGTF = IsoQuant.isoQuantGTF
        File isoQuantOut = IsoQuant.isoQuantOut
        #File isoQuantDenovoGTF = IsoQuantReferenceFree.isoQuantDenovoGTF
        #File isoQuantDenovoOut = IsoQuantReferenceFree.isoQuantDenovoOut
        File stringTieGTF = StringTie.stringTieGTF
        #File stringTieDenovoGTF = StringTieReferenceFree.stringTieDenovoGTF
        File bambuGTF = Bambu.bambuGTF
        File bambuGTFCounts = Bambu.bambuGTFCounts
        File bambuOut = Bambu.bambuOut
        File flairGTF = Flair.flairGTF
        File talonGTF = Talon.talonGTF
        File isoSeqGFF = IsoSeq.isoSeqGFF
        File tamaBED = Tama.tamaBED
        File flamesGFF = Flames.flamesGFF
        File cupcakeGFF = Cupcake.cupcakeGFF
        #File? denovoAnnotationGFFCompareOut = DenovoAnnotationGFFCompare.gffCompareOutput
        File reducedGffCompareOutIsoQuant = ReducedAnnotationGFFCompareIsoQuant.gffCompareOutput
        File reducedGffCompareOutStringTie = ReducedAnnotationGFFCompareStringTie.gffCompareOutput
        File reducedGffCompareOutBambu = ReducedAnnotationGFFCompareBambu.gffCompareOutput
        File reducedGffCompareOutFlair = ReducedAnnotationGFFCompareFlair.gffCompareOutput
        File reducedGffCompareOutTalon = ReducedAnnotationGFFCompareTalon.gffCompareOutput
#        File reducedGffCompareOutIsoSeq = ReducedAnnotationGFFCompareIsoSeq.gffCompareOutput
        File reducedGffCompareOutFlames = ReducedAnnotationGFFCompareFlames.gffCompareOutput
#        File reducedGffCompareOutCupcake = ReducedAnnotationGFFCompareCupcake.gffCompareOutput
        #File? referenceFreeGFFCompareOut = ReferenceFreeGFFCompare.gffCompareOutput
        File isoQuantMonitoringLog = IsoQuant.monitoringLog
        #File isoQuantReferenceFreeMonitoringLog = IsoQuantReferenceFree.monitoringLog
        File stringTieMonitoringLog = StringTie.monitoringLog
        #File stringTieReferenceFreeMonitoringLog = StringTieReferenceFree.monitoringLog
        File bambuMonitoringLog = Bambu.monitoringLog
        File flairMonitoringLog = Flair.monitoringLog
        File talonMonitoringLog = Talon.monitoringLog
        File isoSeqMonitoringLog = IsoSeq.monitoringLog
        File tamaMonitoringLog = Tama.monitoringLog
        File flamesMonitoringLog = Flames.monitoringLog
        File cupcakeMonitoringLog = Cupcake.monitoringLog
        #File denovoAnnotationGFFCompareMonitoringLog = DenovoAnnotationGFFCompare.monitoringLog
        #File referenceFreeGFFCompareMonitoringLog = ReferenceFreeGFFCompare.monitoringLog
        File reducedAnalysisSummary = ReducedAnalysisSummarize.reducedAnalysisSummary
        File reducedAnalysisAccuracyPlots = ReducedAnalysisSummarize.reducedAnalysisAccuracyPlots
    }
}