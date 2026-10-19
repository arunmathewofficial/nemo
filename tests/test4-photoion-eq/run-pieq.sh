#!/bin/bash
#
# Assumes:
# (1) That NEMO has been built in nemo/build/ 
# (2) That the file "data/hii40_test_output.txt" exists for comparison

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

###################################################################
# Logic for whether to test MPI, OpenMP, or hybrid parallelisation
# Default (also with no argument) is to use MPI with np=4
nt=0    # number of threads
mpi=""  # MPI command

if [[ $1 -eq 1 ]]; then
  nt=4
  mpi=""
else
  nt=1
  mpi=""
fi
###################################################################

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"
#export OMP_NUM_THREADS=8
#echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
${mpi} ./photoionization_eq || exit 1

REF_FILE=${script_dir}/../data/pieq_test_refdata.txt
NEW_FILE=output_states_pieq.txt

# Or you can use it in an if statement directly
if python3 ${script_dir}/../compare_data.py $NEW_FILE $REF_FILE 1.0e-4 1.0e-9; then
  echo "Comparison successful - files match"
  echo -e "${GREEN}*** TEST HAS BEEN PASSED ***${NC}"
  exit 0
else
  echo "Comparison failed - files differ or error"
  echo -e "${RED}*** TEST HAS BEEN FAILED ***${NC}"
  exit 1 
fi

