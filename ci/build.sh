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

HERE="$(cd "$(dirname "$0")" && pwd)"
THIS="$(basename "$0")"
NPROC=${NPROC:-2}
WORKD="$(mktemp -d "${THIS}-XXXXXX" -t)"

trap finally ERR
trap finally EXIT

cd "${WORKD}"

export LD_LIBRARY_PATH="${WORKD}/lib:$LD_LIBRARY_PATH"
export PLUGINS_MANIFEST_PATH="${WORKD}/share/plugins"
export ATLAS_DATA_PATH=/usr/local/src/nemo-feedback/atlas-data

ecbuild -S "${HERE}" -DCMAKE_BUILD_TYPE=Debug -DECBUILD_2_COMPAT="ON" -DMPI_ARGS="--oversubscribe"
make -j "${NPROC}"

env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    ctest -j "${NPROC}" -V --output-on-failure --test-dir './nemo-feedback'

exit
