source("tabsEnum.R")

typeSelectionChoices <- reactive({
  if(is.null(full_table_summary())) {
    return(NULL)
  }
  choices <- na.omit(unique(full_table_summary()$Type))
})

setSelectionChoices <- reactive({
  if(is.null(full_table_summary_Legend())) {
    return(NULL)
  }
  choices <- unique(full_table_summary_Legend()$LegendName)
})

stratifierChoices <- reactive({
  if(is.null(full_table_summary_Legend())) {
    return(NULL)
  }
  choices <- unique(full_table_summary_Legend()$Stratifier)
})


engineSelectionChoices <- reactive({
  if(is.null(full_table_summary())) {
    return(NULL)
  }
  choices <- unique(full_table_summary()$Comparison_Engine)
})

output$typeSelectionControl <- renderUI({
  if(!is.null(typeSelectionChoices())) {
    checkboxGroupInput("typeSelect",label=h3("Type"),choices=typeSelectionChoices(),selected=typeSelectionChoices())
  }
})

output$setSelectionControl <- renderUI({
    if(!is.null(setSelectionChoices())) {
      checkboxGroupInput("setSelect",label=h3("Evaluation Sets"),choices=setSelectionChoices(),
                         selected=setSelectionChoices())
    }
})


output$engineSelectionControl <- renderUI({
    if(!is.null(engineSelectionChoices())) {
      checkboxGroupInput("engineSelect",label=h3("Comparison Engine"),choices=engineSelectionChoices(),
                         selected=engineSelectionChoices())
    }
})


output$stratifierSelectionControl <- renderUI({
  if(!is.null(stratifierChoices())) {
    checkboxGroupInput("stratifierSelect",label=h3("Interval Lists"),choices=stratifierChoices(),
                       selected=stratifierChoices())
  }
})


output$downloadButtonControl <- renderUI({
  if(!is.null(full_table_summary())) {
    downloadButton("downloadReport","Download Report")
  }
})
