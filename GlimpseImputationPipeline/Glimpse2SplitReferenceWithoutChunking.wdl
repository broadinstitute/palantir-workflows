version 1.0

workflow Glimpse2SplitReference {
    input {
        # There are two variables required to set define the contigs. The first one, contig_regions, defines the coordinates for the separate chunks,
        # usually this would be an array of ["chr1", "chr2", ..., "chr22", "chrX:1-2781479", "chrX:2781480-155701382", "chrX:155701383-156030895"].
        # Notice the split of chrX into PAR1, non-PAR, PAR2. When constructing the reference panel chunks, however, we need to determine the appropriate
        # reference panel file names. We could do that by parsing the substring before the colon (which WDL doesn't support), or we can just pass
        # another array contig_names_in_reference_panel = ["chr1", "chr2", ..., "chr22", "chrX", "chrX", "chrX"]. Note that contig_regions and
        # contig_names_in_reference_panel must have the same length.
        Array[String] contig_regions
        Array[String] contig_names_in_reference_panel
        Array[String] contig_names_in_genetic_maps
        Array[File] contig_reference_chunks

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

        Int? seed
        Boolean keep_monomorphic_ref_sites = true

        Int preemptible = 1

        # New docker same as old but with bcftools v1.21 instead of v1.16
        String docker = "us.gcr.io/broad-dsde-methods/updated_glimpse_docker:v1.0"
    }

    scatter (i_contig in range(length(contig_regions))) {
        String contig_region = contig_regions[i_contig]
        String contig_name_in_reference_panel = contig_names_in_reference_panel[i_contig]
        String contig_name_in_genetic_maps = contig_names_in_genetic_maps[i_contig]
        File i_reference_chunks = contig_reference_chunks[i_contig]

        String reference_filename = reference_panel_prefix + contig_name_in_reference_panel + reference_panel_suffix
        String genetic_map_filename = genetic_map_path_prefix + contig_name_in_genetic_maps + genetic_map_path_suffix

        call GlimpseSplitReferenceTask {
            input:
                reference_panel = reference_filename,
                reference_panel_index = reference_filename + reference_panel_index_suffix,
                contig = contig_region,
                i_contig = i_contig,
                genetic_map = genetic_map_filename,
                reference_chunks = i_reference_chunks,
                seed = seed,
                keep_monomorphic_ref_sites = keep_monomorphic_ref_sites,
                preemptible = preemptible,
                docker = docker
        }
    }

    output {
        Array[File] reference_chunks = flatten(GlimpseSplitReferenceTask.split_reference_chunks)
    }
}

task GlimpseSplitReferenceTask {
    input {
        String contig
        Int i_contig
        File reference_panel
        File reference_panel_index
        File genetic_map
        File reference_chunks

        Int? seed
        Boolean keep_monomorphic_ref_sites

        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(3.5 * size(reference_panel, "GiB") + size(genetic_map, "GiB") + 100)
        Int preemptible = 0
        String docker
    }

    String reference_output_dir = "reference_output_dir"

    String keep_monomorphic_ref_sites_string = if keep_monomorphic_ref_sites then "--keep-monomorphic-ref-sites" else ""
    command <<<
        set -xeuo pipefail

        NPROC=$(nproc)
        echo "nproc reported ${NPROC} CPUs, using that number as the threads argument for GLIMPSE."

        # Print chunk index to variable
        CONTIGINDEX=$(printf "%04d" ~{i_contig})

        mkdir -p ~{reference_output_dir}

        I_CHUNK=0
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            # Extract coordinates from chunks.txt file
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)

            # Print chunk index to variable
            CHUNKINDEX=$(printf "%04d" $I_CHUNK)

            # Update AC,AN,AF before sending to GLIMPSE to avoid error if these are not up-to-date
            # gcloud storage cat ~{reference_panel} | bcftools +fill-tags /dev/stdin -o updated_ref.vcf.gz -Wtbi -- -t AC,AN,AF
            # export GCS_OAUTH_TOKEN=$(gcloud auth application-default)
            # bcftools view -o current_chunk.vcf.gz -Wtbi --threads $NPROC -r $IRG ~{reference_panel}
            /bin/GLIMPSE2_split_reference --threads ${NPROC} --reference ~{reference_panel} --map ~{genetic_map} --input-region ${IRG} --output-region ${ORG} --output ~{reference_output_dir}/reference_panel_contigindex_${CONTIGINDEX}_chunkindex_${CHUNKINDEX} ~{keep_monomorphic_ref_sites_string} ~{"--seed "+seed}

            # Increase i (and make sure the exit code is zero)
            (( I_CHUNK++ )) || true
        done < ~{reference_chunks}
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        Array[File] split_reference_chunks = glob(reference_output_dir + "/*.bin")
    }
}
