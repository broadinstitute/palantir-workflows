version 1.0

workflow ArrayForDataModel {
  input{
    File file_1
    File file_2
  }
   Array[File] array_of_files = [file_1, file_2]

   output {
    Array[File] array_output = array_of_files
   }
}
