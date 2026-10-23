#!/bin/bash
#
# Assumes:
# (1) That NEMO has been built in nemo/build/ 
# (2) That the file "data/hii40_test_output.txt" exists for comparison

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"

NCPUS=$(getconf _NPROCESSORS_ONLN)
export OMP_NUM_THREADS=$(( NCPUS < 8 ? NCPUS : 8 ))
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
./hii40-lexington #|| exit 1

REF_FILE=${script_dir}/../data/hii40_test_output.txt
NEW_FILE=output_states_hii40.txt

# Or you can use it in an if statement directly
if python3 ${script_dir}/../compare_data.py $NEW_FILE $REF_FILE 3.0e-5 1.0e-10; then
  echo "Comparison successful - files match"
  echo -e "${GREEN}*** TEST HAS BEEN PASSED ***${NC}"
  exit 0
else
  echo "Comparison failed - files differ or error"
  echo -e "${RED}*** TEST HAS BEEN FAILED ***${NC}"
  exit 1 
fi

