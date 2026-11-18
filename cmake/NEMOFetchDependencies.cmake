# Download and build NEMO's dependencies as part of the NEMO build
# (NEMO_FETCH_DEPS=ON), for self-contained builds such as the Python wheels.
#
# - SUNDIALS: CVODE only, serial N_Vector, static libraries, no MPI.
# - fmt and spdlog: static libraries (spdlog uses the external fmt, because
#   NEMO includes <fmt/ranges.h> directly).
# - Boost.Math: header-only, in standalone mode (no other Boost libraries).
#
# Everything is built as position-independent code, so that libnemo.a can be
# linked into a shared library (e.g. a Python extension module).  The
# dependencies are EXCLUDE_FROM_ALL, so nothing of theirs is installed.
#
# To build without network access, point FetchContent at local copies, e.g.
#   -DFETCHCONTENT_SOURCE_DIR_SUNDIALS=/path/to/sundials-7.9.0
#
# Sets: NEMO_SUNDIALS_TARGETS, NEMO_BOOST_TARGET, SUNDIALS_VERSION_MAJOR,
#       SUNDIALS_VERSION_STRING, spdlog_VERSION

cmake_minimum_required(VERSION 3.28)  # FetchContent EXCLUDE_FROM_ALL
include(FetchContent)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(NEMO_SUNDIALS_VERSION "7.9.0")
set(NEMO_FMT_VERSION "12.1.0")
set(NEMO_SPDLOG_VERSION "1.17.0")
set(NEMO_BOOST_MATH_VERSION "1.89.0")

# Settings for the dependencies' own CMake projects.  set(... CACHE ...
# FORCE) would override a user's -D option, so only set normal variables
# (policy CMP0077 makes option() respect them).
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(BUILD_SHARED_LIBS OFF)

# SUNDIALS (option names of SUNDIALS >= 7.9) -------------------------------
set(BUILD_STATIC_LIBS ON)
set(SUNDIALS_ENABLE_ARKODE OFF)
set(SUNDIALS_ENABLE_CVODE ON)
set(SUNDIALS_ENABLE_CVODES OFF)
set(SUNDIALS_ENABLE_IDA OFF)
set(SUNDIALS_ENABLE_IDAS OFF)
set(SUNDIALS_ENABLE_KINSOL OFF)
set(SUNDIALS_ENABLE_FORTRAN OFF)
set(SUNDIALS_ENABLE_MPI OFF)
set(SUNDIALS_ENABLE_OPENMP OFF)
set(SUNDIALS_ENABLE_PTHREAD OFF)
set(SUNDIALS_ENABLE_LAPACK OFF)
set(SUNDIALS_ENABLE_C_EXAMPLES OFF)
set(SUNDIALS_ENABLE_CXX_EXAMPLES OFF)
set(SUNDIALS_ENABLE_EXAMPLES_INSTALL OFF)
set(SUNDIALS_ENABLE_BENCHMARKS OFF)
set(SUNDIALS_ENABLE_PYTHON OFF)

FetchContent_Declare(
  sundials
  URL https://github.com/LLNL/sundials/releases/download/v${NEMO_SUNDIALS_VERSION}/sundials-${NEMO_SUNDIALS_VERSION}.tar.gz
  EXCLUDE_FROM_ALL
)

# fmt and spdlog --------------------------------------------------------------
set(FMT_INSTALL OFF)
set(FMT_TEST OFF)
set(FMT_DOC OFF)
FetchContent_Declare(
  fmt
  URL https://github.com/fmtlib/fmt/archive/refs/tags/${NEMO_FMT_VERSION}.tar.gz
  EXCLUDE_FROM_ALL
)

set(SPDLOG_FMT_EXTERNAL ON)
set(SPDLOG_BUILD_SHARED OFF)
set(SPDLOG_INSTALL OFF)
set(SPDLOG_BUILD_EXAMPLE OFF)
set(SPDLOG_BUILD_TESTS OFF)
FetchContent_Declare(
  spdlog
  URL https://github.com/gabime/spdlog/archive/refs/tags/v${NEMO_SPDLOG_VERSION}.tar.gz
  EXCLUDE_FROM_ALL
)

# Boost.Math (header-only, standalone; not a CMake subproject) --------------
FetchContent_Declare(
  boost_math
  URL https://github.com/boostorg/math/archive/refs/tags/boost-${NEMO_BOOST_MATH_VERSION}.tar.gz
  SOURCE_SUBDIR do-not-add-as-subproject
  EXCLUDE_FROM_ALL
)

FetchContent_MakeAvailable(sundials fmt spdlog boost_math)

# Results ---------------------------------------------------------------------
set(SUNDIALS_VERSION_STRING ${NEMO_SUNDIALS_VERSION})
string(REGEX MATCH "^[0-9]+" SUNDIALS_VERSION_MAJOR ${NEMO_SUNDIALS_VERSION})
set(NEMO_SUNDIALS_TARGETS
    SUNDIALS::cvode_static SUNDIALS::nvecserial_static SUNDIALS::core_static)
foreach(t ${NEMO_SUNDIALS_TARGETS})
  if(NOT TARGET ${t})
    message(FATAL_ERROR "NEMO_FETCH_DEPS: SUNDIALS target ${t} not found")
  endif()
endforeach()

set(spdlog_VERSION ${NEMO_SPDLOG_VERSION})

add_library(nemo_boost_math INTERFACE)
target_include_directories(
  nemo_boost_math SYSTEM INTERFACE ${boost_math_SOURCE_DIR}/include)
target_compile_definitions(nemo_boost_math INTERFACE BOOST_MATH_STANDALONE)
set(NEMO_BOOST_TARGET nemo_boost_math)
