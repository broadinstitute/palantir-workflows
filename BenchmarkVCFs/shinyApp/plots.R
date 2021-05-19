

output$allVariantsPlot <- renderPlot({
  if(!is.null(selectionAllVariantsPlot()) && nrow(selectionAllVariantsPlot())>0) {
      ggplot(selectionAllVariantsPlot(),aes(x=Comparison_Engine,y=value,color=LegendName,shape=Stratifier)) +
        facet_grid(variable~Type) +
        geom_point() +
        scale_y_continuous(limits=range(),expand=c(0,0)) +
        theme(axis.text.x=element_text(angle = 60,hjust=1),
              axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              legend.title=element_blank()
              ) 
  }
},
height=800
)

output$allVariantsHist <- renderPlot({
  if(!is.null(selectionHist()) && nrow(selectionHist())>0) {
      ggplot(selectionHist(),aes(x=LegendName,y=value)) +
        geom_bar(stat="identity",aes(fill=variable)) +
        facet_grid(Type~Comparison_Engine) +
        theme(axis.text.x=element_text(angle = 60,hjust=1),
            axis.text=element_text(size=14),
            axis.title=element_text(size=14),
            legend.title=element_blank()
            ) +
        ylab("counts")
  }
},
height=800
)

output$indelFinePlot <- renderPlot({
  if(!is.null(selectionIndelFine()) && nrow(selectionIndelFine())>0) {
      ggplot(selectionIndelFine(),aes(x=IndelLength,y=value,color=LegendName)) +
        facet_grid(variable~Comparison_Engine) +
        geom_line() +
        geom_point(size=0.7) +
        geom_errorbar(aes(ymin=minus,ymax=plus),width=1.5) +
        scale_y_continuous(limits=range(),expand=c(0,0)) +
        theme(axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              legend.title=element_blank()
              )
  }
},
height=800)

output$indelCoarsePlot <- renderPlot({
  if(!is.null(selectionIndelCoarse()) && nrow(selectionIndelCoarse())>0) {
      ggplot(selectionIndelCoarse(),aes(x=IndelLength,y=value,color=LegendName)) +
        facet_grid(variable~Comparison_Engine) +
        geom_line() +
        geom_point(size=0.7) +
        geom_errorbar(aes(ymin=minus,ymax=plus),width=1.5) +
        scale_y_continuous(limits=range(),expand=c(0,0)) +
        theme(axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              legend.title=element_blank()
              )
  }
},
height=800)