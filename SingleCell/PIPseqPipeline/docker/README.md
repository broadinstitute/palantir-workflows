# Docker and Container Management for Nextflow Pipeline

## Overview

Similar to WDL's `runtime.docker` attribute, Nextflow allows you to specify Docker containers at the **process level** via the `container` directive.

This pipeline uses **two** container images:
- **`qc_container`** (built from `docker/qc/Dockerfile` in this directory): used by `GENERATE_REPORT_DATA`, `CONCATENATE`, `GUIDE_ASSIGNMENT`, and `GENERATE_SUPERSAMPLE_QC`. Bundles pandas/scanpy/anndata/matplotlib plus CRISPAT (installed via git clone + pip, not conda) for guide assignment.
- **`dragen_container`**: used by `DRAGEN_SCRNA`. This is an Illumina-provided DRAGEN image, not built from anything in this repository.

Both are **required** pipeline parameters with no default — you must pass `--qc_container` and `--dragen_container` explicitly.

## Directory Structure

```
docker/
├── README.md              # This file
├── SETUP.md               # Build/push guide for the qc image
├── build_and_push.sh      # Script to build and push the qc image to ECR
└── qc/
    ├── Dockerfile         # qc_container image definition
    └── requirements.txt   # Additional pip-only Python dependencies
```

## Quick Comparison: WDL vs Nextflow

**WDL:**
```wdl
task my_task {
    runtime {
        docker: "my-image:latest"
    }
}
```

**Nextflow:**
```groovy
process MY_PROCESS {
    container "my-image:latest"

    script:
    """
    # your code
    """
}
```

In this pipeline, the container directives reference params rather than hardcoded strings, e.g. `container "${params.qc_container}"` — see `modules/*.nf`.

## Building and Pushing the qc Image

```bash
cd docker
export ECR_REGISTRY=<account-id>.dkr.ecr.<region>.amazonaws.com
export AWS_REGION=<region>
./build_and_push.sh latest
```

See `SETUP.md` for full details, prerequisites, and troubleshooting.

You are not limited to ECR — any registry your executor can pull from works; just build `qc/Dockerfile` and push it with your preferred tool.

## Running

```bash
nextflow run main.nf \
  --qc_container <your-registry>/pipseq-qc:latest \
  --dragen_container <your DRAGEN image> \
  --num_input_cells 10000 \
  --fastq_list fastq_list.csv \
  --supersample_id "Sample_A" \
  --supersample_basename "sample_a" \
  --min_valid_guides 1 \
  --max_valid_guides 2 \
  --ref_tar reference.tar \
  --annotation_file annotation.gtf \
  --outdir results
```

Note that running processes with the `container` directive requires an executor/environment that actually launches containers (e.g. Docker or Singularity enabled in your Nextflow config, or a Kubernetes/ICA executor). This repo does not ship a profile that enables container execution locally — add `docker.enabled = true` (or the Singularity/Podman equivalent) to your own config if you need that.

## Deploying to ICA

DRAGEN in this pipeline is scheduled via Kubernetes pod annotations targeting ICA's FPGA presets (see `pod annotation:` lines in `modules/dragen_scrna.nf`), so the intended deployment target is Illumina Connected Analytics. `export_pipeline_to_ica.py` imports the current git commit of this pipeline into an ICA project as a git-backed Nextflow pipeline; `nextflow_schema.json` drives the parameter form ICA renders from that import, so keep it in sync with any param changes in `main.nf`.

## Updating the Dockerfile

**Conda packages** (preferred for scientific packages):
```dockerfile
# Edit docker/qc/Dockerfile
RUN conda install -y -c conda-forge -c bioconda \
    your-package-name \
    && conda clean -a -y
```

**Pip packages**: edit `docker/qc/requirements.txt` and rebuild.

**crispat itself** is installed via `git clone` + `pip install .` in the Dockerfile (not conda/pip requirements), since it's pulled from a specific fork — see the `git clone` step near the end of `docker/qc/Dockerfile`.

## Testing Locally

```bash
docker build -t my-qc-image:latest -f docker/qc/Dockerfile .
docker run -it my-qc-image:latest python -c "import pandas, scanpy, crispat; print('OK')"
```

## Troubleshooting

**Problem: "Container not found"**
- Ensure the image is pushed to a registry your executor can pull from
- Check image name/tag spelling and that `--qc_container`/`--dragen_container` were actually passed

**Problem: "Permission denied" for scripts**
- Scripts in `bin/` are automatically made available in `PATH` by Nextflow — no need to `chmod +x` inside the container

**Problem: "Module not found" in Python**
- Add the missing package to `docker/qc/Dockerfile` (or `requirements.txt`) and rebuild/push the `qc_container` image

**Problem: Different behavior locally vs. in the cloud**
- Check platform compatibility: add `--platform linux/amd64` to `docker build` if building on Apple Silicon
- Verify your registry region/credentials match what your executor expects

## Additional Resources

- Nextflow containers: https://www.nextflow.io/docs/latest/container.html
- AWS ECR: https://docs.aws.amazon.com/AmazonECR/latest/userguide/what-is-ecr.html
- CRISPAT Documentation: https://github.com/pinellolab/CRISPAT
