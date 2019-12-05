#!/usr/bin/env Rscript

library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("Please supply 2 arguments, input monitoring.log file and output pdf file.")
}
lines=readLines(args[1])
pdfout=args[2]

theme_update(text = element_text(family="Helvetica"), axis.text=element_text(family="Helvetica"),legend.text=element_text(family="Helvetica"))

prevDate=NA

process_monitor_line = function(x) {
  if(grepl("^\\[.*\\]$",x) ) {
    prevDate<<-as.POSIXct(x, format="[%a %b %d %H:%M:%S UTC %Y]")
    return (NULL)
  }
  if(grepl("CPU usage",x)) {
    text=gsub("\\* (.*): (.*)%",'data.frame(variable="\\1",value1="\\2",unit1="%",value2="",unit2="")',x)
    
  } else {
  text=gsub("\\* (.*): (.*?) (.*?) (\\d+(\\.\\d+)?)(.*?)$",'data.frame(variable="\\1",value1="\\2",unit1="\\3",value2="\\4",unit2="\\6")',x)
  }
  vec=eval(parse(text=text))
  ret = cbind(vec,data.frame(date=prevDate))
}

headline=grep("--- Runtime Information ---",lines,fixed = T)[1]
lines=lines[(headline+1):length(lines)]
data=ldply(lines,process_monitor_line)
recasted=melt(id.vars="date",dcast(melt(data = data,id.vars = c("date","variable"),variable.name = "variable2"),formula=date~variable+variable2))

vars=gsub("_value.","",unique(subset(recasted,grepl("value",variable),select="variable")$variable))
units=transform(ddply(subset(recasted,grepl("unit",variable)),.(variable),function(x) c(value=unique(x$value))),variable=gsub("_unit.","",variable),index=gsub(".*_unit","",variable))
values_only=subset(recasted,grepl("value",variable))%>%transform(variable=gsub("_value.","",variable),index=gsub(".*_value","",variable))

pdf(pdfout,width = 8,height = 12)
ggplot(values_only)+geom_line(aes(x=date,y=as.numeric(as.character(value))))+
  facet_wrap(~variable+index,scale="free_y")+scale_x_datetime()+ 
  geom_text(data=units,x=-Inf,y=Inf, hjust=0,vjust=1,color="red",aes(label=value),family="Helvetica")+labs(y="value")
dev.off()



