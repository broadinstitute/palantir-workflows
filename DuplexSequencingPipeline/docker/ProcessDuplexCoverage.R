library(plyr)
library(ggplot2)
#options(stringsAsFactors=F)
readDuplexCoverage <- function(coverageFile) {
  read.table(coverageFile, header=T, sep = "\t")
}

processResults <- function(coverageFile) {
   duplexCoverage = readDuplexCoverage(coverageFile)
  
   results = ddply(duplexCoverage, names(duplexCoverage), .progress = "text", function(x)  {
      s = strsplit(as.character(x[, names(duplexCoverage)[5]]), " ")[[1]]
  
      A = as.numeric(gsub(".:", "", s[1]))
      C = as.numeric(gsub(".:", "", s[2]))
      G = as.numeric(gsub(".:", "", s[3]))
      T = as.numeric(gsub(".:", "", s[4]))
      N = as.numeric(gsub(".:", "", s[5]))
  
      ref = max(A, T, C, G)
      bases = c(A, T, C, G)
      alt = max(bases[-which.max(bases)])
      data.frame(A, T, C, G, N, ref, alt)
   })
   results
   #results[results$alt > 0, c("Locus", "A", "T", "C", "G", "N", "ref", "alt")]
}

generateDepthFigures <- function(baseName, rawDepthFile, startStopDepthFile, duplexDepthFile) {
  rawDepth = processResults(rawDepthFile)
  startStopDepth = processResults(startStopDepthFile)
  duplexDepth = processResults(duplexDepthFile)

  # Create plot of start-stop deduplicated depth by duplex depth
  ggplot(merge(duplexDepth, startStopDepth, by="Locus"), aes(x = Total_Depth.x, y = Total_Depth.y)) + xlab("Duplex Depth") + ylab("Start-Stop Depth") + geom_point()
  ggsave(paste0(baseName, "-StartStopVsDuplex.pdf"))
  
  # Create plot of raw read depth vs duplex depth
  ggplot(merge(duplexDepth, rawDepth, by="Locus"), aes(x = Total_Depth.x, y = Total_Depth.y)) + xlab("Duplex Depth") + ylab("Raw Depth") + geom_point()
  ggsave(paste0(baseName, "-RawVsDuplex.pdf"))
  
  # Create plot of raw read depth vs start-stop depth
  ggplot(merge(rawDepth,startStopDepth, by="Locus"), aes(x = Total_Depth.x, y = Total_Depth.y)) + xlab("Raw Depth") + ylab("Start-Stop Depth") + geom_point()
  ggsave(paste0(baseName, "-StartStopVsRawDepth.pdf"))
  
  meanRawDepth = mean(rawDepth$ref + rawDepth$alt)
  meanStartStopDepth = mean(startStopDepth$ref + startStopDepth$alt)
  meanDuplexDepth = mean(duplexDepth$ref + duplexDepth$alt)
  meanDepth = data.frame(meanRawDepth, meanStartStopDepth, meanDuplexDepth)
  write.table(meanDepth, paste0(baseName, ".depth.txt"), sep = "\t", row.names = F, quote = FALSE)  
}
