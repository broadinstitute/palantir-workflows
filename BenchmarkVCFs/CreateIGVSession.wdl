version 1.0

# creates an IGV session
# given a list of IGV compatible file paths
workflow CreateIGVSession{
    input {
        Array[String] input_bams
        Array[String] input_vcfs
        String reference_version
        String output_name
    }

    call WriteXMLfile {
        input:
            input_files= flatten([input_bams, input_vcfs]),
            reference_version=reference_version,
            file_name=output_name
    }
}

task WriteXMLfile {
    input {
        Array[String] input_files
        String reference_version
        String file_name

        Array[String]? input_names=""
        Array[String] input_names_prefix = if defined(input_names) then prefix('-n ', select_first([input_names])) else []
    }
    command {
        bash /usr/writeIGV.sh ~{reference_version} ~{sep=" " input_files} ~{sep=" " input_names_prefix}  > "~{file_name}.xml"
    }
    runtime {
        docker: "quay.io/mduran/generate-igv-session_2:v1.0"
    }
    output {
        File igv_session = "${file_name}.xml"
    }
}