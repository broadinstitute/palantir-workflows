# PALANTIR-WORKFLOWS

[![CircleCI](https://circleci.com/gh/broadinstitute/palantir-workflows.svg?style=svg)](https://circleci.com/gh/broadinstitute/palantir-workflows)

Utility workflows used by the DSP's Palantir team.  This repository should be used to manage frequently used utility workflows for the team, and facilitate their use on [Terra](https://app.terra.bio/) through [Dockstore](https://dockstore.org/).

**Remember, this is a public repository, so anything you put in this repo is publicly viewable.**


## Testing Workflows

All workflows should have associated tests.
In order to add tests, you should add a test workflow to the `test` directory. 
The test workflow should call the workflow you are testing, and (preferably) compare the outputs to those expected. 
Input JSONs for the test workflow must be placed in a directory whose name is the same as the test workflow, with `.wdl` replaced by `_json`.
So, the test directory structure will be built like this:

```bash
+-- palantir-workflows
|   +-- test
|   |   +-- MyWorkflow
|   |   |   +-- my_test_workflow.wdl
|   |   |   +-- my_test_workflow_json
|   |   |   |   +-- test_input_1.json
|   |   |   |   +-- test_input_2.json
+++++++++++++++++
```
