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

    # Shard the VCF file into chunks and process for GLIMPSE
    scatter (interval in read_lines(contig_reference_chunks)) {
        call ShardVcf {
            input:
                vcf = reference_filename,
                vcf_index = reference_filename_index,
                interval=interval,

        }

        call GlimpseSplitReferenceTask {
            input:
                reference_panel = ShardVcf.vcf_chunk,
                reference_panel_index = ShardVcf.vcf_chunk_index,
                contig = contig_name,
                interval = interval,
                genetic_map = genetic_map_filename,
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

        String interval
        Boolean add_allele_info = true

        Int disk_size = ceil(2.5 * size(vcf, "GiB") + 100)
    }

    command <<<
        set -xueo pipefail

        INTERVAL=$(echo "~{interval}" | awk '{ print $3 }')
        # Use bcftools to split the VCF file into chunks
        if [[ "~{add_allele_info}" == "true" ]]; then
            bcftools view -t $INTERVAL ~{vcf} --threads $(nproc) | bcftools +fill-tags - -o "chunk.vcf.gz" -Wtbi -- -t AC,AN,AF
        else
            bcftools view -t $INTERVAL ~{vcf} --threads $(nproc) -o "chunk.vcf.gz" -Wtbi
        fi
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/updated_glimpse_docker:v1.0"
        memory: "8 GB"
        cpu: 4
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

    output {
        File vcf_chunk = "chunk.vcf.gz"
        File vcf_chunk_index = "chunk.vcf.gz.tbi"
    }
}

task GlimpseSplitReferenceTask {
    input {
        String contig
        File reference_panel
        File reference_panel_index
        File genetic_map
        String interval

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

        # Print contig index to variable
        CONTIGINDEX="~{contig}"

        # Make chunk index from interval
        INTERVAL=$(echo "~{interval}" | awk '{ print $3 }')
        CHUNKINDEX=$(echo "${INTERVAL}" | tr ":" "-")

        mkdir -p reference_output_dir

        /bin/GLIMPSE2_split_reference \
            --threads ${NPROC} \
            --reference ~{reference_panel} \
            --map ~{genetic_map} \
            --input-region ${IRG} \
            --output-region ${ORG} \
            --output reference_output_dir/reference_panel_contigindex_${CONTIGINDEX}_chunkindex_${CHUNKINDEX} \
            ~{keep_monomorphic_ref_sites_string} \
            ~{"--seed " + seed}
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
