#!/usr/bin/env bash

set -e 

test_file=${1}
test_cdl="${test_file::-3}.cdl"
ref_cdl=${2}

ncdump ${test_file} >> ${test_cdl}

diff ${test_cdl} ${ref_cdl}
