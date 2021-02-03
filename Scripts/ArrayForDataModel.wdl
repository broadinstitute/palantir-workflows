workflow fof_usage_wf {
   File file_of_files
   Array[File] array_of_files = read_lines(file_of_files)

   output {
    Array[File] array_output = array_of_files
   }
}
