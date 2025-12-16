version 1.0

workflow testOOM {
    call testOOM_t
}

task testOOM_t {
    command <<<
        echo "Killed"
        exit 1
    >>>

    runtime {
        docker:"ubuntu:24.04"
        memory: "10 GB"
        maxRetries: 2
    }
}