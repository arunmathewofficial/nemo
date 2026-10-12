#!/bin/sh

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"
build_dir="${script_dir}/build"

mkdir -p ${build_dir}
pushd ${build_dir}

cmake \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}" \
    "${script_dir}"

make -j 4
popd


