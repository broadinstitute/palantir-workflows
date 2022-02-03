# Build docker image from Dockerfile in this directory
# This is to help keep track of versions, etc.
# Wraps around ImputationPipeline/build_push_docker.sh.

dockerfile_directory=tidyverse_kableextra_docker
image_version=v1.0.0 # as of Jan 13 2022

wd=$(pwd)
cd "$(dirname $0)" || exit

../build_push_docker.sh \
  --directory ${dockerfile_directory} \
  --image-version-tag ${image_version} \
#  --dry-run


cd "${wd}" || exit