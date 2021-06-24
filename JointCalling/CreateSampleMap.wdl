version 1.0

workflow CreateSampleMapFile{
    
    input {
        Array[String] sample_names
        Array[File] gvcfs
    }


    
    call CreateSampleMap {
        input:
            sample_names=sample_names,
            input_gvcfs=gvcfs
    }
    output {
        File sample_map_file = CreateSampleMap.output_map
    }
}


task CreateSampleMap {
    input {
        Array[String] sample_names
        Array[File] input_gvcfs
    }

    parameter_meta {
        input_gvcfs: {
          localization_optional: true
        }
    }

    command <<<

        set -xe

        python << CODE
        gvcfs = ['~{sep="','" input_gvcfs}']
        samples = ['~{sep="','" sample_names}']

        if len(gvcfs) != len(samples):
          print(f"error! len(gvcfs)={len(gvcfs)} which is different from len(sample_names)={len(samples)}. Quiting!")    
          exit(1)

        with open("inputs.list", "w") as fi:
          for i in range(len(gvcfs)):
            fi.write(samples[i] + "\t" + gvcfs[i] + "\n") 

        CODE

    >>>
    runtime {
        docker: "python:3.6"
        memory: "7 GB"
        cpu: "2"
        disks: "local-disk " + 50 + " HDD"
    }
    output {
        File output_map = "inputs.list"
    }
}