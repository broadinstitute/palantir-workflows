# ICA Compatibility Setup

This pipeline has been configured to be compatible with the [`nextflow-to-icav2-config`](https://github.com/keng404/nextflow-to-icav2-config) tool for deployment to Illumina Connected Analytics (ICA).

## Changes Made

The following modifications were made to enable ICA compatibility:

### 1. Created `nextflow_schema.json`
- Defines all pipeline parameters with types, descriptions, and validation rules
- Used by the ICA tool to generate the parameters XML file
- Organizes parameters into logical groups (input/output, pipeline options, etc.)

### 2. Created `conf/base.ica.config`
- Contains ICA-specific pod annotations for compute instance selection
- Defines default container for all processes
- Sets process-specific CPU/memory requirements
- Configures Docker settings for ICA

### 3. Updated `nextflow.config`
- Added manifest closure with pipeline metadata
- Included reference to `conf/base.ica.config`
- Added ICA-specific parameters (`ica_smoke_test`)
- Enabled parameter validation

### 4. Added Container Directives
- All processes now have explicit container directives
- Container: `us.gcr.io/broad-dsde-methods/pipseq-qc:latest`
- Ensures ICA knows which Docker image to use for each process

## Using the ICA Conversion Tool

### Prerequisites

1. Install the tool's Docker image:
```bash
docker pull keng404/nextflow-to-icav2-config:0.0.7
```

2. Obtain an ICA API key (see [ICA documentation](https://help.connected.illumina.com/account-management/platform-home#manage-existing-api-keys))

3. Create a project in ICA (via GUI or CLI)

### Running the Conversion

From the `palantir-workflows/SingleCell` directory:

```bash
docker run -itv $(pwd):$(pwd) --platform=linux/amd64 keng404/nextflow-to-icav2-config:0.0.7 /bin/bash
```

Inside the container:

```bash
cd $(pwd)/PIPseqPipeline

# Run the conversion
Rscript /path/to/nf-core.conversion_wrapper.R \
  --pipeline-dirs $(pwd) \
  --staging-directory $(pwd) \
  --run-scripts /path/to/nextflow-to-icav2-config \
  --intermediate-copy-template /path/to/nextflow-to-icav2-config/dummy_template.txt \
  --create-pipeline-in-ica \
  --api-key-file /path/to/api-key.txt \
  --ica-project-name "YOUR_PROJECT_NAME" \
  --nf-core-mode \
  --in-docker
```

Or use individual steps:

#### Step 1: Generate parameters XML
```bash
Rscript create_xml/nf-core.json_to_params_xml.R --json nextflow_schema.json
```

#### Step 2: Update configs for ICA
```bash
Rscript ica_nextflow_config.test.R --config-file nextflow.config
```

#### Step 3: Add ICA modifications to Nextflow scripts
```bash
Rscript develop_mode.downstream.R \
  --config-file nextflow.config \
  --nf-script main.nf \
  --other-workflow-scripts modules/generate_report_data.nf \
  --other-workflow-scripts modules/concatenate.nf \
  --other-workflow-scripts modules/guide_assignment.nf \
  --other-workflow-scripts modules/generate_supersample_qc.nf
```

#### Step 4: Create pipeline in ICA
```bash
Rscript nf-core.create_ica_pipeline.R \
  --nextflow-script main.nf \
  --workflow-language nextflow \
  --parameters-xml parameters.xml \
  --nf-core-mode \
  --ica-project-name "YOUR_PROJECT_NAME" \
  --pipeline-name "pipseq-qc" \
  --api-key-file /path/to/api-key.txt
```

## What the Tool Will Do

The `nextflow-to-icav2-config` tool will:

1. **Parse your configuration files** and update them with ICA-specific settings
2. **Add workflow.onError handlers** for better troubleshooting in ICA
3. **Modify bin/ script references** so ICA can properly execute them
4. **Generate parameters XML** from `nextflow_schema.json` for ICA's input form
5. **Add additional ICA-specific edits** for smooth execution

## ICA Compute Instances

Pod annotations in `conf/base.ica.config` map to these instance types:

| Pod Annotation Value | CPUs | Memory |
|---------------------|------|--------|
| `standard-small` | 2 | 4 GB |
| `standard-medium` | 4 | 8 GB |
| `standard-large` | 8 | 16 GB |
| `standard-xlarge` | 16 | 32 GB |

See [ICA documentation](https://help.ica.illumina.com/project/p-flow/f-pipelines#definition) for full list.

## Container Requirements

Before deploying to ICA, ensure your Docker container is:
1. Pushed to a registry accessible by ICA (GCR, ECR, Docker Hub)
2. Contains all required dependencies (Python, pandas, scanpy, anndata, crispat, etc.)
3. Has AWS CLI v2 installed (if using S3 storage)

Update the container path in `conf/base.ica.config` if using a different registry.

## Testing

After conversion, test locally:
```bash
Rscript testing_pipelines/test_nextflow_script.R \
  --nextflow-script main.nf \
  --docker-image nextflow/nextflow:22.04.3 \
  --nextflow-config nextflow.config
```

## Notes

- ICA currently supports Nextflow versions `22.04.3` and `20.10.0`
- The `ica_smoke_test` parameter can be set to `true` to run test configurations
- All processes must have container directives for ICA execution
- Pod annotations are required for ICA to allocate compute resources
