version 1.0

task PerformPCA {
  input {
    File bim
    File bed
    File fam
    String basename
    Int n_pcs = 20
    Int mem = 8
    Int nthreads = 16
  }

  Int mem_gb = mem * 1024 / 2

  # again, based on Wallace commands
  command <<<
    cp ~{bim} ~{basename}.bim
    cp ~{bed} ~{basename}.bed
    cp ~{fam} ~{basename}.fam

    ~/flashpca/flashpca \
      --bfile ~{basename} \
      -n ~{nthreads} \
      -d ~{n_pcs} \
      --memory ~{mem_gb} \
      --outpc ~{basename}.pc \
      --outpve ~{basename}.pc.variance \
      --outload ~{basename}.pc.loadings \
      --outmeansd ~{basename}.pc.meansd
  >>>

  output {
    File pcs = "~{basename}.pc"
    File pc_variance = "~{basename}.pc.variance"
    File pc_loadings = "~{basename}.pc.loadings"
    File mean_sd = "~{basename}.pc.meansd"
    File eigenvectors = "eigenvectors.txt"
    File eigenvalues = "eigenvalues.txt"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
    disks: "local-disk 400 HDD"
    memory: mem + " GB"
  }
}

# This projects the array dataset using the previously generated PCs, using flashPCA
task ProjectArray {
  input {
    File bim
    File bed
    File fam
    File pc_loadings
    File pc_meansd
    String basename
    Int mem = 8
    Int nthreads = 16
  }

  Int mem_gb = mem * 1024 / 2

  command <<<
    cp ~{bim} ~{basename}.bim
    cp ~{bed} ~{basename}.bed
    cp ~{fam} ~{basename}.fam

    cp ~{pc_loadings} loadings.txt
    cp ~{pc_meansd} meansd.txt

    # Check if .bim file, pc loadings, and pc meansd files have the same IDs
    # 1. extract IDs, removing first column of .bim file and first rows of the pc files
    awk '{print $2}' ~{basename}.bim > bim_ids.txt
    awk '{print $1}' loadings.txt | tail -n +2 > pcloadings_ids.txt
    awk '{print $1}' meansd.txt | tail -n +2 > meansd_ids.txt

    diff bim_ids.txt pcloadings_ids.txt > diff1.txt
    diff bim_ids.txt meansd_ids.txt > diff2.txt
    diff pcloadings_ids.txt meansd_ids.txt > diff3.txt

    if [[ -s diff3.txt ]]
    then
    echo "PC loadings file and PC means file do not contain the same IDs; check your input files and run again."
    exit 1
    fi

    # check if diff files are not empty
    if [[ -s diff1.txt || -s diff2.txt ]]
    then
    echo "IDs in .bim file are not the same as the IDs in the PCA files; check that you have the right files and run again."
    exit 1
    fi

    ~/flashpca/flashpca \
      --bfile ~{basename} \
      --numthreads ~{nthreads} \
      --memory ~{mem_gb} \
      --project \
      --inmeansd meansd.txt \
      --outproj projections.txt \
      --inload loadings.txt \
      -v
  >>>

  output {
    File projections = "projections.txt"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
    disks: "local-disk 400 HDD"
    memory: mem + " GB"
  }
}

task ArrayVcfToPlinkDataset {
  input {
    File vcf
    File pruning_sites
    File? subset_to_sites
    String basename
    Int mem = 8
  }

  Int disk_space =  3 * ceil(size(vcf, "GB")) + 20

  command <<<
    /plink2 \
      --vcf ~{vcf} \
      --extract-intersect ~{pruning_sites} ~{subset_to_sites} \
      --allow-extra-chr \
      --set-all-var-ids @:#:\$1:\$2 \
      --new-id-max-allele-len 1000 missing \
      --out ~{basename} \
      --make-bed \
      --rm-dup force-first
  >>>

  output {
    File bed = "~{basename}.bed"
    File bim = "~{basename}.bim"
    File fam = "~{basename}.fam"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    disks: "local-disk " + disk_space + " HDD"
    memory: mem + " GB"
  }
}