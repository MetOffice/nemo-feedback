#!/usr/bin/env bash

set -euo pipefail

test_file="$1"
ref_cdl="$2"

diff -u <(ncdump "${test_file}") "${ref_cdl}"
