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

# latest JCSDA images use spack to install the environment
#. /etc/profile.d/z10_spack_environment.sh

if [ -z ${LD_LIBRARY_PATH:-} ]; then
  export LD_LIBRARY_PATH="${WORKD}/lib"
else
  export LD_LIBRARY_PATH="${WORKD}/lib:${LD_LIBRARY_PATH}"
fi
#export PLUGINS_MANIFEST_PATH="${WORKD}/share/plugins"

rm -f "${HERE}/nemo-feedback"
ln -s '..' "${HERE}/nemo-feedback"
ecbuild -S "${HERE}" -DCMAKE_BUILD_TYPE=Debug -DMPI_ARGS="--oversubscribe"
make -j "${NPROC}"

env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
   ctest -j "${NPROC}" -V --output-on-failure --test-dir "./nemo-feedback"

exit
