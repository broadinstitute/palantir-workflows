version 1.0

workflow CombineTables {
    input {
        Array[File] tables
    }

    call CombineTablesScript {
        input:
            tables=tables
    }

    output {
        File combined_table = CombineTablesScript.combined_table
    }
}

task CombineTablesScript {
    input {
        Array[File] tables
        String comment_char = "#"
        String output_name = "combined_table"
    }

    command <<<
        set -xueo

        python << CODE
        import pandas as pd

        combined_df = pd.DataFrame()
        for table in ["~{sep="\", \"" tables}"]:
            combined_df = pd.concat([combined_df, pd.read_csv(table, sep="\t", comment="~{comment_char}")])

        combined_df.to_csv("~{output_name}.tsv", sep="\t", index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk " + ceil(2.5 * size(tables, "GB") + 20) + " HDD"
        cpu: 2
        memory: "4GB"
    }

    output {
        File combined_table = "~{output_name}.tsv"
    }
}