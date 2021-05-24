
observeEvent(input$table_rows_selected, {
  rows <- input$table_rows_selected
  names <- lapply(rows,function(x) dirname(gsub("gs://","https://console.cloud.google.com/storage/browser/",table_summary_Legend()$IGV_Session[x])))
  browseToIGVSessionName(names)
  }
)

observeEvent(input$tableFine_rows_selected, {
  rows <- input$tableFine_rows_selected
  names <- lapply(rows,function(x) dirname(gsub("gs://","https://console.cloud.google.com/storage/browser/",table_summary_fine_Legend()$IGV_Session[x])))
  browseToIGVSessionName(names)
  }
)

observeEvent(input$tableCoarse_rows_selected, {
  rows <- input$tableCoarse_rows_selected
  names <- lapply(rows,function(x) dirname(gsub("gs://","https://console.cloud.google.com/storage/browser/",table_summary_coarse_Legend()$IGV_Session[x])))
  browseToIGVSessionName(names)
  }
)

observeEvent(input$tableInsertion_Deletion_rows_selected, {
  rows <- input$tableInsertion_Deletion_rows_selected
  names <- lapply(rows,function(x) dirname(gsub("gs://","https://console.cloud.google.com/storage/browser/",table_summary_insertion_deletion()$IGV_Session[x])))
  browseToIGVSessionName(names)
  }
)

browseToIGVSessionName <- function(names) {
  for (name in names) {
    js$browseURL(name)
  }
}