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

rm -f "${HERE}/nemo-feedback"
ln -s '..' "${HERE}/nemo-feedback"

ecbuild -S "${HERE}" -DCMAKE_BUILD_TYPE=Debug -DECBUILD_2_COMPAT="ON" -DMPI_ARGS="--oversubscribe"
make -j "${NPROC}"
env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    ctest -j "${NPROC}" --output-on-failure --test-dir './nemo-feedback'

exit
