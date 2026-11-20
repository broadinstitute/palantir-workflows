# PALANTIR-WORKFLOWS

[![Test Status](https://github.com/broadinstitute/palantir-workflows/actions/workflows/run_tests.yaml/badge.svg?branch=main)](https://github.com/broadinstitute/palantir-workflows/actions/workflows/run_tests.yaml/badge.svg?branch=main) [![Dockerfiles Status](https://github.com/broadinstitute/palantir-workflows/actions/workflows/test_dockerfiles.yaml/badge.svg?branch=main)](https://github.com/broadinstitute/palantir-workflows/actions/workflows/test_dockerfiles.yaml/badge.svg?branch=main)

Utility workflows used by the DSP's Palantir team.  This repository should be used to manage frequently used utility workflows for the team, and facilitate their use on [Terra](https://app.terra.bio/) through [Dockstore](https://dockstore.org/).

**Remember, this is a public repository, so anything you put in this repo is publicly viewable.**


## Testing Workflows

Automated WDL testing is implemented using [watt](https://github.com/rickymagner/watt).
To add tests, update  [test/watt_config.yml](test/watt_config.yml).
See the watt documentation for usage details.  Tests are run using github actions, controlled by [.github/workflows/run_tests.yaml](.github/workflows/run_tests.yaml).  Automated testing also validates all WDLs found in the repository, regardless of whether tests have been added for the particular WDL, using `womtool validate`.  Automated testing is performed on all PR's, as well as any pushes to the `main` branch.

During testing, WDLs are run on GCP Batch in the `broad-gp-hydrogen` project, with the execution bucket `gs://broad-gp-hydrogen-palantir-workflows-test-execution`.  Generally, input and expected output files are stored in `gs://broad-gp-hydrogen-palantir-workflows-test-data`.  The stdout from watt will be visible in the github actions UI, and cromwell log files can be found in a zipped artifact named `cromwell_logs`, also in the github actions UI.   

### CI cost

Test VMs and disks are labeled through Cromwell's `google_labels` workflow option (written by the `Write Cromwell options` step in [run_tests.yaml](.github/workflows/run_tests.yaml) and passed to watt with `--cromwell-options`): `palantir-ci=true`, `ci-trigger=<event, e.g. pull-request>` and `ci-run=run-<GitHub run id>`. The test buckets should also carry `palantir-ci=true` (`gcloud storage buckets update gs://<bucket> --update-labels=palantir-ci=true`) so storage cost appears under the same label. Cloud Billing data is exported to the BigQuery dataset `broad-gp-hydrogen:billing_export`; queries for daily cost, cost per CI run, and cost per Cromwell workflow are in [test/billing](test/billing).

## Using the Dockstore Github App to Automatically Update Workflows in Dockstore/Terra
Workflows registered in Dockstore can be automatically synced when changes are pushed to this repo by adding their information to `.dockstore.yml`. 
In this way, a change pushed to a branch in this repo can be automatically propagated into any Terra workspaces using the workflow. 
Details can be found at Dockstore's Github App [documentation](https://docs.dockstore.org/en/develop/getting-started/github-apps/github-apps.html).
