#!/usr/bin/env bash

set -euo pipefail

test_file="$1"
ref_cdl="$2"

diff -u "${ref_cdl}" <(ncdump "${test_file}")
