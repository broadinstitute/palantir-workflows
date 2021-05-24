#Compute +/- 1-sigma binomial errors using Wilson Score interval (https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval)

GetWilsonUpIndel <- function(variable,value,TP_Eval,TP_Base,FP,FN) {
  if (variable=="F1_Score") {
    return(value)
  }
  if (variable=="Precision") {
    n=TP_Eval+FP
    return((value+1/(2*n))/(1+1/n)+1/(1+1/n)*sqrt(value*(1-value)/n+1/(4*n^2)))
  }
  if (variable=="Recall"){
    n=TP_Base+FN
    return((value+1/(2*n))/(1+1/n)+1/(1+1/n)*sqrt(value*(1-value)/n+1/(4*n^2)))
  }
}
GetWilsonDownIndel <- function(variable,value,TP_Eval,TP_Base,FP,FN) {
  if (variable=="F1_Score") {
    return(value)
  }
  if (variable=="Precision") {
    n=TP_Eval+FP
    return((value+1/(2*n))/(1+1/n)-1/(1+1/n)*sqrt(value*(1-value)/n+1/(4*n^2)))
  }
  if (variable=="Recall"){
    n=TP_Base+FN
    return((value+1/(2*n))/(1+1/n)-1/(1+1/n)*sqrt(value*(1-value)/n+1/(4*n^2)))
  }
}
library(reshape2)

full_table_summary <- reactive({
  if(is.null(input$inputFiles)) {
    return(NULL)
  }
  ret <-read.csv(input$inputFiles$datapath)
  ret$Stratifier[is.na(ret$Stratifier)]="N/A"
  ret
})

full_table_summary_Legend <- reactive({
  if(is.null(full_table_summary())) {
    return(NULL)
  }
  table <- full_table_summary()
  table$LegendName<-paste(as.character(table$Name),"vs",as.character(table$Truth_Set))
  return(table)
})

table_summary_insertion_deletion <- reactive({
  if(is.null(full_table_summary())) {
    return(NULL)
  }
  
  ret <-subset(full_table_summary(),Summary_Type%in%c("insertion","deletion"))
  
  if(nrow(ret)>0) {
    return(ret)
  }
  return(NULL)
  
})

table_summary_fine_Legend <- reactive({
  if(is.null(full_table_summary_Legend())) {
    return(NULL)
  }
  ret <-subset(full_table_summary_Legend(),Summary_Type=="indel_fine")
  if(nrow(ret)>0) {
    return(ret)
  }
  return(NULL)
})



table_summary_Legend <- reactive({
  if(is.null(full_table_summary_Legend())) {
    return(NULL)
  }
  ret <-subset(full_table_summary_Legend(),Summary_Type=="summary")
  if(nrow(ret)>0) {
    return(ret)
  }
  return(NULL)
})

baseVariables=c("Name", "Truth_Set", "Comparison_Engine","Type","LegendName","Stratifier","IGV_Session")
extraVariables=c("IndelLength","TP_Eval","TP_Base","FP","FN")
performanceVariables=c("Recall","Precision","F1_Score")

melted_table <- reactive({
  table<- table_summary_Legend()
  if(is.null(table)) {
    return(NULL)
  }
  ret<-melt(table,id.vars=baseVariables,measure.vars=performanceVariables)
})

melted_table_hist <- reactive({
  if(is.null(table_summary_Legend())) {
    return(NULL)
  }
  ret <-melt(table_summary_Legend(),id.vars=baseVariables,measure.vars=c("TP_Eval","FN","FP","UNK"))
})

melted_table_fine <- reactive({
  if(is.null(table_summary_fine_Legend())) {
    return(NULL)
  }
  ret<-melt(table_summary_fine_Legend(),id.vars=c(baseVariables,extraVariables) ,measure.vars=performanceVariables)
  ret$plus<-mapply(GetWilsonUpIndel,as.character(ret$variable),ret$value,ret$TP_Eval,ret$TP_Base,ret$FP,ret$FN)
  ret$minus<-mapply(GetWilsonDownIndel,as.character(ret$variable),ret$value,ret$TP_Eval,ret$TP_Base,ret$FP,ret$FN)
  return(ret)
})

table_summary_coarse_Legend <- reactive({
  if(is.null(full_table_summary_Legend())) {
    return(NULL)
  }
  ret <-subset(full_table_summary_Legend(),Summary_Type=="indel_coarse")
  if(nrow(ret)>0) {
    return(ret)
  }
  return(NULL)
})

melted_table_coarse <- reactive({
  if(is.null(table_summary_coarse_Legend())) {
    return(NULL)
  }
  ret<-melt(table_summary_coarse_Legend(),id.vars=c(baseVariables,extraVariables),measure.vars=performanceVariables)
  ret$plus<-mapply(GetWilsonUpIndel,as.character(ret$variable),ret$value,ret$TP_Eval,ret$TP_Base,ret$FP,ret$FN)
  ret$minus<-mapply(GetWilsonDownIndel,as.character(ret$variable),ret$value,ret$TP_Eval,ret$TP_Base,ret$FP,ret$FN)
  return(ret)
})

