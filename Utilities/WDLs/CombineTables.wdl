version 1.0

workflow CombineTables {
    input {
        Array[File] tables
        String output_name = "combined_table"
        Array[String] extra_column_names = []    # Names of extra columns to add to output
        Array[String] extra_column_values = []    # Values to put in extra columns; one value per column name
    }

    call CombineTablesScript {
        input:
            tables=tables,
            extra_column_names=extra_column_names,
            extra_column_values=extra_column_values,
            output_name=output_name
    }

    output {
        File combined_table = CombineTablesScript.combined_table
    }
}

task CombineTablesScript {
    input {
        Array[File] tables
        Array[String] extra_column_names
        Array[String] extra_column_values

        String comment_char = "#"
        String output_name
    }

    command <<<
        set -xueo

        python << CODE
        import pandas as pd

        combined_df = pd.DataFrame()
        for table in ["~{sep="\", \"" tables}"]:
            combined_df = pd.concat([combined_df, pd.read_csv(table, sep="\t", comment="~{comment_char}")])

        name_val_dict = dict(zip(["~{sep="\", \"" extra_column_names}"], ["~{sep="\", \"" extra_column_values}"]))
        for key, value in name_val_dict.items():
            combined_df[key] = value

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