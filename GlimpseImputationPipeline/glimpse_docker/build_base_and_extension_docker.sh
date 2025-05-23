#!/usr/bin/env bash

# This branch is using the following configuration.
# To replicate this branch, run this script with the following arguments:
# build_base_and_extension_docker.sh -r "https://github.com/kachulis/GLIMPSE.git" -b ck_checkpoint_clean -t "glimpse_build_test"

set -Eeuo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)

usage() {
  cat <<EOF
Usage: $(basename "${BASH_SOURCE[0]}") [-h] [-p] -t docker_tag -r glimpse_repo -b glimpse_branch

Build docker image to be used in the Glimpse2Imputation WDL. This script will clone the
specified GLIMPSE repo on the specified branch and build the docker image. Subsequently,
it will build the extension docker image for the WDL based on the previously built
GLIMPSE docker image, tag it, and optionally push the image.

Available options:

-h, --help      Print this help and exit.
-p, --push      If set, push image.
-y, --yes       If set, don't ask for confirmation.
-t, --tag       Tag for the newly created docker image. If the tag does not contain a colon, it will be determined automatically based on the repo and commit hash.
-r, --repo      Repo URL for the GLIMPSE2 base image (e.g. https://github.com/odelaneau/GLIMPSE.git).
-b, --branch    Branch name for the GLIMPSE2 base image (e.g. master).
EOF
  exit
}

cleanup() {
  trap - SIGINT SIGTERM ERR EXIT
  rm -rf "${script_dir}/glimpse_base"
}

setup_colors() {
  if [[ -t 2 ]] && [[ -z "${NO_COLOR-}" ]] && [[ "${TERM-}" != "dumb" ]]; then
    NOFORMAT='\033[0m' RED='\033[0;31m' GREEN='\033[0;32m' ORANGE='\033[0;33m' BLUE='\033[0;34m' PURPLE='\033[0;35m' CYAN='\033[0;36m' YELLOW='\033[1;33m'
  else
    NOFORMAT='' RED='' GREEN='' ORANGE='' BLUE='' PURPLE='' CYAN='' YELLOW=''
  fi
}

msg() {
  echo >&2 -e "${1-}"
}

end_program() {
  local msg=$1
  local code=${2-1} # default exit status 1
  msg "$msg"
  msg "Run $(basename "${BASH_SOURCE[0]}") --help for usage information."
  exit "$code"
}

parse_params() {
  # default values of variables set from params
  push=0
  yes=0

  while :; do
    case "${1-}" in
    -h | --help) usage ;;
    --no-color) NO_COLOR=1 ;;
    -p | --push) push=1 ;;
    -y | --yes) yes=1 ;;
    -t | --tag)
      tag="${2-}"
      shift
      ;;
    -r | --repo)
      repo="${2-}"
      shift
      ;;
    -b | --branch)
      branch="${2-}"
      shift
      ;;
    -?*) end_program "Unknown option: $1" ;;
    *) break ;;
    esac
    shift
  done

  # check required params and arguments
  [[ -z "${tag-}" ]] && end_program "Missing required parameter: tag"
  [[ -z "${repo-}" ]] && end_program "Missing required parameter: repo"
  [[ -z "${branch-}" ]] && end_program "Missing required parameter: branch"

  return 0
}

parse_params "$@"
setup_colors

pushing=$(if [ "${push}" -eq "1" ]; then echo "pushing"; else echo "not pushing"; fi)

tag_string=$(if [[ "$tag" != *":"* ]]; then echo "${YELLOW}automatically${NOFORMAT} tagging it based on the repo name and commit hash"; else echo "${YELLOW}manually${NOFORMAT} tagging it as ${YELLOW}${tag}${NOFORMAT}"; fi)

msg "Building docker image based on GLIMPSE repo ${YELLOW}${repo}${NOFORMAT} on branch ${YELLOW}${branch}${NOFORMAT} and ${tag_string} and ${YELLOW}${pushing}${NOFORMAT} it."

if [ "${yes}" -eq "0" ]; then
    read -p "Continue (y/n)? " choice
    case "$choice" in 
        y|Y ) msg "Continuing...";;
        * ) end_program "Aborting.";;
    esac
fi

if [ -d "${script_dir}/glimpse_base" ]; then
    end_program "The glimpse_base subdirectory already exist. Please remove it before building the docker."
fi

trap cleanup SIGINT SIGTERM ERR EXIT

if ! docker images > /dev/null; then
    end_program "Couldn't determine docker images, probably because the docker deamon isn't running. Please start docker and restart this script."
    exit 1
fi

if docker images | grep "temp_glimpse_base" > /dev/null; then
    end_program "temp_glimpse_base Docker image exists already. Please remove it before building to ensure that the correct base image is used: docker image rm temp_glimpse_base"
    exit 1
fi

git clone $repo --branch $branch --single-branch ${script_dir}/glimpse_base

if [[ "$tag" != *":"* ]]; then
    prefix="https://github.com/"
    suffix="/GLIMPSE.git"
    if [[ "$repo" != "$prefix"*"$suffix" ]]; then
        end_program "repo does not match the expected format. When providing only a Docker path and not a tag, the URL must follow this format: $prefix<GITHUB_NAMESPACE>$suffix"
        exit 1
    fi

    github_namespace=${repo#"$prefix"}
    github_namespace=${github_namespace%"$suffix"}
    
    glimpse_commit_hash=$(git -C "${script_dir}/glimpse_base" rev-parse --short HEAD)
    tag="${tag}:${github_namespace}_${glimpse_commit_hash}"
fi

docker build --platform linux/amd64 -t temp_glimpse_base ${script_dir}/glimpse_base
docker build --platform linux/amd64 -t ${tag} ${script_dir}

if [ "$push" -eq "1" ]; then
    msg "Pushing docker image to ${tag}"
    docker push ${tag}
fi
