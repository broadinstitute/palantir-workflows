# PALANTIR-WORKFLOWS

[![CircleCI](https://circleci.com/gh/broadinstitute/palantir-workflows.svg?style=svg)](https://circleci.com/gh/broadinstitute/palantir-workflows)

Utility workflows used by the DSP's Palantir team.  This repository should be used to manage frequently used utility workflows for the team, and facilitate their use on [Terra](https://app.terra.bio/) through [Dockstore](https://dockstore.org/).

**Remember, this is a public repository, so anything you put in this repo is publicly viewable.**


## Testing Workflows

All (WDL) workflows should have associated tests.
Tests can be added by adding jobs to the `validate-and-test` (circleci) workflow in `.circleci/config.yml`.
Each (WDL) workflow should be at least validated with womtool by adding the following job to the `validate-and-test` (circleci) workflow.
	
	- validate:
		wdl: "path/to/workflow/myWorkflow.wdl"
		name: "validate-myWorkflow"

If possible, each (WDL) workflow should also be tested by running some plumbing test data through to make sure that it succeeds.
This can be done by additionally adding the following job to the `validate-and-test` (circleci) workflow

	- test:
		wdl: "path/to/workflow/myWorkflow.wdl"
		input-json: "path/to/workflow/myWorkflow.wdl.test.json"
		name: "test-myWorkflow"
		requires:
			-validate-myWorkflow

Input paths in the input json should be absolute, with the consideration that the `palantir-workflows` directory will map to `/home/circleci/project` in the test environment.
Tests can also consist of test (WDL) workflow, which import the (WDL) workflow to be tested, run it, and compare the outputs to expected outputs.
See `test/BenchmarkVCFs/test_BenchmarkVCFs.wdl` for example.