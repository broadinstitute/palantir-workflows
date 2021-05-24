library(shinyjs)

# define js function for opening urls in new tab/window
js_code <- "
shinyjs.browseURL = function(url) {
window.open(url,'_blank');
}
"

source("tabsEnum.R")

ui <- fluidPage(
  useShinyjs(),
  extendShinyjs(text = js_code, functions = 'browseURL'),

  sidebarPanel(
    fileInput("inputFiles","Data Files",multiple=TRUE,accept = c(".csv")),
    uiOutput("downloadButtonControl"),
    conditionalPanel(
      condition = paste0("input.tabSelected!=",toString(tabs$Histogram)," & input.tabSelected<",toString(tabs$Table)),
      sliderInput("range", label = h3("Y-Range"), min = 0, 
                  max = 1, value = c(0, 1))
    ),
    conditionalPanel(
      condition = paste0("input.tabSelected<",tabs$Indels_fine_plot),
      uiOutput("typeSelectionControl")
    ),
    conditionalPanel(
      condition = paste0("input.tabSelected!=",toString(tabs$Histogram)," & input.tabSelected<",toString(tabs$Table)),
      checkboxGroupInput("metricSelect",label=h3("Metric"),choices=list("Recall"="Recall","Precision"="Precision","F1 Score"="F1_Score"),
                         selected=c("Recall","Precision","F1_Score"))
    ),
    conditionalPanel(
      condition = paste0("input.tabSelected<",toString(tabs$Table)),
      uiOutput("engineSelectionControl")
    ),
    conditionalPanel(
      condition = paste0("input.tabSelected<",toString(tabs$Table)),
      uiOutput("setSelectionControl")
    ),
    
    conditionalPanel(
      condition = paste0("input.tabSelected<",toString(tabs$Table)),
      uiOutput("stratifierSelectionControl")
    )
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Metrics",value=tabs$Metrics,plotOutput(outputId = "allVariantsPlot",click="plot_click_allVariantsPlot")),
      tabPanel("Histogram",value=tabs$Histogram,plotOutput(outputId = "allVariantsHist")),
      tabPanel("Indels (Fine Binning)",value=tabs$Indels_fine_plot,plotOutput(outputId = "indelFinePlot",click="plot_click_finePlot")),
      tabPanel("Indels (Coarse Binning)",value=tabs$Indels_coarse_plot,plotOutput(outputId = "indelCoarsePlot",click="plot_click_coarsePlot")),
      tabPanel("Table",value=tabs$Table,DT::dataTableOutput(outputId = "table")),
      tabPanel("Indel Length Table (Fine Binning)",value=tabs$Table_fine,DT::dataTableOutput(outputId = "tableFine")),
      tabPanel("Indel Length Table (Coarse Binning)",value=tabs$Table_coarse,DT::dataTableOutput(outputId = "tableCoarse")),
      tabPanel("Insertion/Deletion Table",value=tabs$Table_insert_delete,DT::dataTableOutput(outputId = "tableInsertion_Deletion")),
      id="tabSelected"
    )
  )
)