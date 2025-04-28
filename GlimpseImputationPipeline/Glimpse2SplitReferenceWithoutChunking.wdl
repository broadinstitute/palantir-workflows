version 1.0

workflow Glimpse2SplitReference {
    input {
        # There are two variables required to set define the contigs. The first one, contig_regions, defines the coordinates for the separate chunks,
        # usually this would be an array of ["chr1", "chr2", ..., "chr22", "chrX:1-2781479", "chrX:2781480-155701382", "chrX:155701383-156030895"].
        # Notice the split of chrX into PAR1, non-PAR, PAR2. When constructing the reference panel chunks, however, we need to determine the appropriate
        # reference panel file names. We could do that by parsing the substring before the colon (which WDL doesn't support), or we can just pass
        # another array contig_names_in_reference_panel = ["chr1", "chr2", ..., "chr22", "chrX", "chrX", "chrX"]. Note that contig_regions and
        # contig_names_in_reference_panel must have the same length.
        String contig_name
        File contig_reference_chunks

        # Example path for chr1: gs://bucket/to/panel/reference_panel_chr1_merged.bcf(.csi)
        # reference_panel_prefix = "gs://bucket/to/panel/reference_panel_"
        # reference_panel_suffix = "_merged.bcf"
        # reference_panel_index_suffix = ".csi"
        String reference_panel_prefix
        String reference_panel_suffix
        String reference_panel_index_suffix

        # Same format as reference_panel_pre-/suffix
        String genetic_map_path_prefix
        String genetic_map_path_suffix

        Int seed = 42
        Boolean keep_monomorphic_ref_sites = true

        Int lines_per_chunk = 5
        Boolean add_allele_info = true

        # New docker same as old but with bcftools v1.21 instead of v1.16
        String docker = "us.gcr.io/broad-dsde-methods/updated_glimpse_docker:v1.0"
    }

    String reference_filename = reference_panel_prefix + contig_name + reference_panel_suffix
    String reference_filename_index = reference_filename + reference_panel_index_suffix
    String genetic_map_filename = genetic_map_path_prefix + contig_name + genetic_map_path_suffix

    call ShardVcf {
        input:
            vcf = reference_filename,
            vcf_index = reference_filename_index,
            ref_chunks = contig_reference_chunks,
            lines_per_chunk = lines_per_chunk,
            contig_name = contig_name
    }

    # Keep the three: vcf_chunk, vcf_chunk_index, ref_chunk together with a counter at end (chunk.right.right)
    scatter (chunk in zip(zip(ShardVcf.vcf_chunks, ShardVcf.vcf_chunks_indices), zip(ShardVcf.ref_chunks, range(length(ShardVcf.ref_chunks))))) {
        if (add_allele_info) {
            call AddAlleleInfo {
                input:
                    vcf = chunk.left.left,
                    vcf_index = chunk.left.right,
                    shard_number = chunk.right.right
            }
        }

        call GlimpseSplitReferenceTask {
            input:
                reference_panel = select_first([AddAlleleInfo.updated_vcf, chunk.left.left]),
                reference_panel_index = select_first([AddAlleleInfo.updated_vcf_index, chunk.left.right]),
                contig = contig_name,
                genetic_map = genetic_map_filename,
                reference_chunks = chunk.right.left,
                seed = seed,
                keep_monomorphic_ref_sites = keep_monomorphic_ref_sites,
                docker = docker
        }
    }


    output {
        Array[File] reference_chunks = flatten(GlimpseSplitReferenceTask.split_reference_chunks)
    }
}

task ShardVcf {
    input {
        File vcf
        File vcf_index

        File ref_chunks
        Int lines_per_chunk = 5
        Int max_jobs = 10
        String contig_name

        Int disk_size = ceil(2.5 * size(vcf, "GiB") + 100)
    }

    command <<<
        set -xueo pipefail

        python3 << CODE
        import pandas as pd

        df = pd.read_csv("~{ref_chunks}", sep="\t", header=None, names=["index", "contig", "IR", "OR", "cm", "c1", "c2", "c3"])
        df["region"] = df["IR"].apply(lambda x: x.split(":")[1])
        df["start"] = df["region"].apply(lambda x: int(x.split("-")[0]))
        df["end"] = df["region"].apply(lambda x: int(x.split("-")[1]))

        # Split dataframe into chunks of size lines_per_chunk
        chunks = [df[i:i + ~{lines_per_chunk}] for i in range(0, df.shape[0], ~{lines_per_chunk})]

        # Create final intervals using first start and final end values
        final_intervals = [f'~{contig_name}:{d["start"].values[0]}-{d["end"].values[-1]}' for d in chunks]

        # Write final_intervals to file
        with open("shard_intervals.txt", "w") as f:
            f.write("\n".join(final_intervals))

        # Write each chunk to a separate file
        for i, chunk in enumerate(chunks):
            chunk.to_csv(f"chunk_{i}.txt", sep="\t", index=False, header=False)

        CODE

        # Use bcftools to split the VCF file into chunks
        I_CHUNK=0

        # Parallelize the IO
        MAX_JOBS=~{max_jobs}  # Limit concurrent jobs
        JOB_COUNT=0

        while IFS= read -r interval || [[ -n "$interval" ]]; do
            bcftools view -r "$interval" ~{vcf} -Oz -o "chunk_${I_CHUNK}.vcf.gz" -Wtbi --threads $(nproc) &

            # Move lower to start counting chunks at 0
            (( I_CHUNK++ )) || true
            (( JOB_COUNT++ )) || true
            if (( JOB_COUNT >= MAX_JOBS )); then
                wait -n  # Wait for any job to finish
                (( JOB_COUNT-- )) || true
            fi
        done < shard_intervals.txt

        wait  # Wait for all jobs to finish

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/updated_glimpse_docker:v1.0"
        memory: "8 GB"
        cpu: 4
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

    output {
        Array[File] vcf_chunks = glob("chunk_*.vcf.gz")
        Array[File] vcf_chunks_indices = glob("chunk_*.vcf.gz.tbi")
        Array[String] intervals = read_lines("shard_intervals.txt")
        Array[File] ref_chunks = glob("chunk_*.txt")
    }
}

task AddAlleleInfo {
    input {
        File vcf
        File vcf_index

        Int shard_number

        Int disk_size = ceil(2.5 * size(vcf, "GiB") + 100)
    }

    command <<<
        set -xueo pipefail

        # Use bcftools to add allele information to the VCF file
        bcftools +fill-tags ~{vcf} -Oz -o "updated_chunk_~{shard_number}.vcf.gz" -Wtbi -- -t AC,AN,AF
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/updated_glimpse_docker:v1.0"
        memory: "8 GB"
        cpu: 4
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

    output {
        File updated_vcf = "updated_chunk_~{shard_number}.vcf.gz"
        File updated_vcf_index = "updated_chunk_~{shard_number}.vcf.gz.tbi"
    }
}

task GlimpseSplitReferenceTask {
    input {
        String contig
        File reference_panel
        File reference_panel_index
        File genetic_map
        File reference_chunks

        Int? seed
        Boolean keep_monomorphic_ref_sites

        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(2.2 * size(reference_panel, "GiB") + size(genetic_map, "GiB") + 100)
        String docker
    }

    String keep_monomorphic_ref_sites_string = if keep_monomorphic_ref_sites then "--keep-monomorphic-ref-sites" else ""

    command <<<
        set -xeuo pipefail

        NPROC=$(nproc)
        echo "nproc reported ${NPROC} CPUs, using that number as the threads argument for GLIMPSE."

        # Print chunk index to variable
        CONTIGINDEX="~{contig}"

        mkdir -p reference_output_dir

        I_CHUNK=0
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            # Extract coordinates from chunks.txt file
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)

            # Print chunk index to variable
            CHUNKINDEX=$(printf "%04d" $I_CHUNK)

            /bin/GLIMPSE2_split_reference \
                --threads ${NPROC} \
                --reference ~{reference_panel} \
                --map ~{genetic_map} \
                --input-region ${IRG} \
                --output-region ${ORG} \
                --output reference_output_dir/reference_panel_contigindex_${CONTIGINDEX}_chunkindex_${CHUNKINDEX} \
                ~{keep_monomorphic_ref_sites_string} \
                ~{"--seed "+seed}

            # Increase i (and make sure the exit code is zero)
            (( I_CHUNK++ )) || true
        done < ~{reference_chunks}
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: 0
    }

    output {
        Array[File] split_reference_chunks = glob("reference_output_dir/*.bin")
    }
}
