#!/bin/bash
#
# Test NEMO with more than one radiation source: run the HII40 Lexington
# problem with one star ("single"), with the star split into two identical
# stars with half of the luminosity each ("split"), and with a second,
# different star that does not illuminate any cell ("absent").  The three
# results must agree (no reference data is needed).
#
# Assumes that it is run in the directory with the hii40-lexington
# executable.

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"

NCPUS=$(getconf _NPROCESSORS_ONLN)
export OMP_NUM_THREADS=$(( NCPUS < 8 ? NCPUS : 8 ))
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"

# The exit status of hii40-lexington is not used: it is non-zero if the
# comparison with the published benchmark values fails or cannot be made,
# as in run-hii40-test.sh.  Check that the output file was written instead.
for mode in single split absent; do
  rm -f output_2src_${mode}.txt
  ./hii40-lexington ${mode} output_2src_${mode}.txt
  if [ ! -s output_2src_${mode}.txt ]; then
    echo -e "${RED}*** No output from mode ${mode}: TEST HAS BEEN FAILED ***${NC}"
    exit 1
  fi
done

# Relative and absolute tolerances for the comparison with "single".
# "absent" should give the same numbers.  For "split" the sum of the two
# half-rates differs from the single rate by round-off, which the ODE
# solver (relative tolerance 1e-6) amplifies.  Measured differences: below
# 1e-10 inside the ionised region; in the three outermost (neutral) cells,
# which are at the minimum temperature, 1e-8 in temperature and up to
# 1.5e-5 in the trace He1+ fraction; 1e-13 (absolute) in trace ions.
status=0
compare() {
  mode=$1
  if python3 ${script_dir}/../compare_data.py output_2src_${mode}.txt \
      output_2src_single.txt $2 $3; then
    echo "Mode ${mode}: same result as with a single star"
  else
    echo "Mode ${mode}: result differs from that with a single star"
    status=1
  fi
}
compare split 1.0e-4 1.0e-11
compare absent 1.0e-12 1.0e-20

if [ ${status} -eq 0 ]; then
  echo -e "${GREEN}*** TEST HAS BEEN PASSED ***${NC}"
else
  echo -e "${RED}*** TEST HAS BEEN FAILED ***${NC}"
fi
exit ${status}
