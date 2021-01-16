output$downloadReport <- downloadHandler(
  filename="report.pdf",
  content = function(file) {
    tempReport <- file.path(tempdir(),"report.Rmd")
    file.copy("report.Rmd", tempReport, overwrite = TRUE)
    params <- list(range=range(),
                   selectionAllVariants=selectionAllVariantsPlot(),
                   selectionHist=selectionHist(),
                   selectionIndelFine=selectionIndelFine(),
                   selectionIndelCoarse=selectionIndelCoarse()
    )
    
    withProgress(message='Downloading Report', value=0, {
    rmarkdown::render(tempReport, output_file = file,
                      params = params,
                      envir = new.env(parent = globalenv())
    )
    })
  }
)