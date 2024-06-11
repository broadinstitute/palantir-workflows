# PALANTIR-WORKFLOWS

[![CircleCI](https://circleci.com/gh/broadinstitute/palantir-workflows.svg?style=svg)](https://circleci.com/gh/broadinstitute/palantir-workflows)

Utility workflows used by the DSP's Palantir team.  This repository should be used to manage frequently used utility workflows for the team, and facilitate their use on [Terra](https://app.terra.bio/) through [Dockstore](https://dockstore.org/).

**Remember, this is a public repository, so anything you put in this repo is publicly viewable.**


## Testing Workflows

Automated WDL testing is implemented using [watt](https://github.com/rickymagner/watt).
To add tests, update  [test/watt_config.yml](test/watt_config.yml).
See the watt documentation for usage details.
BLA BLA BLA GCP
Automated testing on CircleCI also validates every WDL in the repo using the `validate` tool from `womtool`. 

## Using the Dockstore Github App to Automatically Update Workflows in Dockstore/Terra
Workflows registered in Dockstore can be automatically synced when changes are pushed to this repo by adding their information to `.dockstore.yml`. 
In this way, a change pushed to a branch in this repo can be automatically propagated into any Terra workspaces using the workflow. 
Details can be found at Dockstore's Github App [documentation](https://docs.dockstore.org/en/develop/getting-started/github-apps/github-apps.html).
