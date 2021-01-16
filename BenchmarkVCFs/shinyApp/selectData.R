selectionAllVariantsPlot <- reactive({
  if(is.null(melted_table())) {
    return(NULL) 
  }
  sel<-subset(melted_table(),Type%in%input$typeSelect & variable%in%input$metricSelect & LegendName%in%input$setSelect & Stratifier%in%input$stratifierSelect
              & Comparison_Engine%in%input$engineSelect)
  sel <- as.data.frame(lapply(sel, function (x) if (is.factor(x)) factor(x) else x))
})

selectionHist <- reactive({
  if(is.null(melted_table_hist())) {
    return(NULL)
  }
  sel<-subset(melted_table_hist(),variable%in%c("TP_Eval","FP","FN","UNK") & LegendName%in%input$setSelect & Type%in%input$typeSelect & Stratifier%in%input$stratifierSelect
              & Comparison_Engine%in%input$engineSelect)
})

range <- reactive({
  ret <- c(input$range[1],input$range[2])
})

selectionIndelFine <- reactive({
  if(is.null(melted_table_fine())) {
    return(NULL)
  }
  sel <- subset(melted_table_fine(), variable%in%input$metricSelect & LegendName%in%input$setSelect & Stratifier%in%input$stratifierSelect
                & Comparison_Engine%in%input$engineSelect)
})

selectionIndelCoarse <- reactive({
  if(is.null(melted_table_coarse())) {
    return(NULL)
  }
  sel <- subset(melted_table_coarse(), variable%in%input$metricSelect & LegendName%in%input$setSelect & Stratifier%in%input$stratifierSelect
                & Comparison_Engine%in%input$engineSelect)
})