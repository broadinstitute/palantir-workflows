observeEvent(input$plot_click_allVariantsPlot, {
  rows <- nearPoints(selectionAllVariantsPlot(),input$plot_click_allVariantsPlot)
  browseToIGVSession(rows)
}
)

observeEvent(input$plot_click_finePlot, {
  rows <- nearPoints(selectionIndelFine(),input$plot_click_finePlot)
  browseToIGVSession(rows)
}
)

observeEvent(input$plot_click_coarsePlot, {
  rows <- nearPoints(selectionIndelCoarse(),input$plot_click_coarsePlot)
  browseToIGVSession(rows)
}
)

browseToIGVSession <- function(rows) {
  names <- lapply(rows$IGV_Session,function(x) dirname(gsub("gs://","https://console.cloud.google.com/storage/browser/",x)))
  for (name in names) {
    js$browseURL(name)
  }
} 
