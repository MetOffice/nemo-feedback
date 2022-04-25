#!/bin/bash
set -euo pipefail

finally() {
    trap '' ERR
    trap '' EXIT
    if [[ -d "${WORKD:-}" ]]; then
        cd /
        rm -fr "${WORKD}"
    fi
}

apt-get -y install libhdf5-serial-dev \
apt-get -y install libnetcdf-dev libnetcdff-dev libnetcdf-c++4-dev \
apt-get -y install lz4

HERE="$(cd "$(dirname "$0")" && pwd)"
THIS="$(basename "$0")"
NPROC=${NPROC:-2}
WORKD="$(mktemp -d "${THIS}-XXXXXX" -t)"

trap finally ERR
trap finally EXIT

cd "${WORKD}"

rm -f "${HERE}/nemo-feedback"
ln -s '..' "${HERE}/nemo-feedback"
ecbuild -S "${HERE}" -DMPI_ARGS="--oversubscribe"
make -j "${NPROC}"
env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    ctest -j "${NPROC}" --output-on-failure --test-dir './nemo-feedback'

exit
