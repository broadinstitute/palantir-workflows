# Build docker image from Dockerfile in this directory
# This is to help keep track of versions, etc.
# Wraps around PRS/build_push_docker.sh.

dockerfile_directory=tidyverse_kableextra_docker
image_version=v1.1.0 # as of Feb 8 2022

wd=$(pwd)
cd "$(dirname $0)" || exit

../build_push_docker.sh \
  --directory ${dockerfile_directory} \
  --image-version-tag ${image_version} \
#  --dry-run


cd "${wd}" || exit