identifierOrder <- c("Name","Truth_Set","Comparison_Engine","Stratifier")
metricOrder <- c("Recall","Precision","F1_Score")
countOrder <- c("TP_Base","TP_Eval","FP","FN")

output$table <- DT::renderDataTable({
  if(!is.null(table_summary_Legend())) {
    formatRound(DT::datatable(
                              subset(table_summary_Legend(),select=-c(LegendName,IGV_Session,Summary_Type,IndelLength))[,c(3,4,5,11,12,1,2,6,7,8,9,10,13)],
                              rownames=FALSE,
                              selection="single",
                              options=list(
                                order=list(list(0,"asc"),list(1,"asc"),list(2,"asc"),list(3,"asc")),
                                columnDefs = list(list(className = 'dt-center', targets = c(0,1,2,3)))
                                )
                              ),
                digits=3,columns=c(6,7,8)
                )
    }
})

output$tableFine <- DT::renderDataTable({
  if(!is.null(table_summary_fine_Legend())) {
    formatRound(DT::datatable(
                                subset(table_summary_fine_Legend(),select=-c(LegendName,IGV_Session,Type,UNK,Summary_Type))[,c(3,4,5,11,12,1,2,6,7,8,9,10)],
                                rownames=FALSE,
                                selection="single",
                                options=list(
                                  order=list(list(0,"asc"),list(1,"asc"),list(2,"asc"),list(3,"asc"),list(4,"asc")),
                                  columnDefs = list(list(className = 'dt-center', targets = c(0,1,2,3,4)))
                                  )
                              ),
                digits=3,columns=c(6,7,8)
                )
    }
})

output$tableCoarse <- DT::renderDataTable({
  if(!is.null(table_summary_coarse_Legend())) {
    formatRound(DT::datatable(
                                subset(table_summary_coarse_Legend(),select=-c(LegendName,IGV_Session,Type,UNK,Summary_Type))[,c(3,4,5,11,12,1,2,6,7,8,9,10)],
                                rownames=FALSE,
                                selection="single",
                                options=list(
                                  order=list(list(0,"asc"),list(1,"asc"),list(2,"asc"),list(3,"asc"),list(4,"asc")),
                                  columnDefs = list(list(className = 'dt-center', targets = c(0,1,2,3,4)))
                                  )
                              ),
                                digits=3,columns=c(6,7,8)
                )  
    }
})

output$tableInsertion_Deletion <- DT::renderDataTable({
  if(!is.null(table_summary_insertion_deletion())) {
    formatRound(DT::datatable(
                                subset(table_summary_insertion_deletion(),select=-c(IGV_Session,Type,UNK,IndelLength))[,c(3,4,5,11,12,1,2,6,7,8,9,10)],
                                rownames=FALSE,
                                selection="single",
                                options=list(
                                  order=list(list(0,"asc"),list(1,"asc"),list(2,"asc"),list(3,"asc"),list(4,"asc")),
                                  columnDefs = list(list(className = 'dt-center', targets = c(0,1,2,3,4)))
                                  )
                              ),
                                digits=3,columns=c(6,7,8)
                )  
    }
})