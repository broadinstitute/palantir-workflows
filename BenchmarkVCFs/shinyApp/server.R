library(ggplot2)
library(reshape2)
library(DT)

options(DT.options = list(pageLength=25))

server <- function(input, output) {
  source("loadData.R",local=TRUE)
  source("selectData.R",local=TRUE)
  source("uiControl.R",local=TRUE)
  source("plots.R",local=TRUE)
  source("plotClicks.R",local=TRUE)
  source("tableClicks.R",local=TRUE)
  source("tables.R",local=TRUE)
  source("download.R",local=TRUE)
}