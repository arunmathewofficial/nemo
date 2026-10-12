#!/bin/sh

script="${BASH_SOURCE[0]:-${(%):-%x}}"
script_dir="$( cd "$( dirname "${script}" )" >/dev/null 2>&1 && pwd )"
build_dir="${script_dir}/build"

mkdir -p ${build_dir}
pushd ${build_dir}

CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}"
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}" \
    -DNEMO_OMP=ON \
    -DOpenMP_CXX_LIB_NAMES="libomp" \
    -DOpenMP_libomp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib \
    -DOpenMP_CXX_INCLUDE_DIR=/opt/homebrew/opt/libomp/include/ \
    "${script_dir}"

make -j 4
popd

