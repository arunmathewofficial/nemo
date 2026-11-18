cmake_minimum_required(VERSION 3.15.0)

project(
    NEMO
    VERSION 1.4.0
    LANGUAGES CXX
)


# Sane default install directories
include(GNUInstallDirs)

# User provided install directories
set(
    NEMO_INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR} CACHE STRING
    "The install directory for NEMO executables."
)

set(
    NEMO_INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR} CACHE STRING
    "The install directory for NEMO libraries."
)

set(
    NEMO_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE STRING
    "The install directory for NEMO headers."
)

set(
    NEMO_INSTALL_DATADIR ${CMAKE_INSTALL_DATADIR} CACHE STRING
    "The install directory for NEMO data (eg. cmake scripts)."
)

option(
    NEMO_PYTHON
    "Build the Python module nemo_plasma._core (sets NEMO_FETCH_DEPS=ON, \
NEMO_BUILD_TESTS=OFF, NEMO_OMP=OFF).  Normally set by pip via pyproject.toml."
    OFF
)
if(NEMO_PYTHON)
  set(NEMO_FETCH_DEPS ON)
  set(NEMO_BUILD_TESTS OFF)
  set(NEMO_OMP OFF)
endif()

option(
    NEMO_OMP
    "Build NEMO with OpenMP. Default is disbled."
    OFF
)

option(
    NEMO_FETCH_DEPS
    "Download and build SUNDIALS (CVODE), spdlog, fmt and Boost.Math as part \
of the NEMO build (static, no MPI), instead of using installed packages. \
NEMO is then not installed.  Used for the Python package."
    OFF
)

option(
    NEMO_BUILD_TESTS
    "Build the NEMO tests."
    ON
)

#  This is required on OSX so that the Boost code will compile successfully
set(CMAKE_CXX_STANDARD 14)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

if(NEMO_OMP)
  add_definitions(-DNEMO_OMP)
endif(NEMO_OMP)


# Find dependant libraries
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

if(NEMO_OMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  include_directories(${OpenMP_CXX_INCLUDE_DIR})
endif(NEMO_OMP)

if(NEMO_FETCH_DEPS)
  include(NEMOFetchDependencies)
  set(NEMO_MPI_TARGET "")
else()
  find_package(Boost 1.73.0 REQUIRED)

  find_package(SUNDIALS REQUIRED)
  #message(STATUS "SUNDIALS version = ${SUNDIALS_VERSION_MAJOR}.${SUNDIALS_VERSION_MINOR}")
  #message(STATUS "SUNDIALS include = ${SUNDIALS_INCLUDE_DIR}/../")

  find_package(MPI REQUIRED)  # Homebrew sundials needs open MPI

  find_package(spdlog REQUIRED)
  add_definitions("-DSPDLOG_COMPILED_LIB")
  # Get spdlog version
  message(STATUS "spdlog version: ${spdlog_VERSION}")

  # Create an interface target that abstracts the spdlog+fmt combination
  add_library(spdlog_deps INTERFACE)

  if(spdlog_VERSION VERSION_LESS "1.13.0" OR APPLE)
    message(STATUS "using external fmt b/c old spdlog or apple")
    find_package(fmt REQUIRED)
    message(STATUS "fmt version = ${fmt_VERSION}")
    target_compile_definitions(spdlog_deps INTERFACE SPDLOG_FMT_EXTERNAL)
    target_link_libraries(spdlog_deps INTERFACE spdlog::spdlog fmt::fmt)
    # Suppress deprecation warnings from fmt/spdlog version mismatch in Homebrew
    target_compile_options(spdlog_deps INTERFACE -Wno-deprecated-declarations)
  else()
    target_link_libraries(spdlog_deps INTERFACE spdlog::spdlog)
  endif()

  set(NEMO_BOOST_TARGET Boost::boost)
  set(NEMO_MPI_TARGET MPI::MPI_CXX)
  set(NEMO_SUNDIALS_TARGETS ${SUNDIALS_LIBRARY})
  if(${SUNDIALS_VERSION_STRING} VERSION_GREATER_EQUAL 7.0.0)
    list(APPEND NEMO_SUNDIALS_TARGETS ${SUNDIALS_sundials_core_LIBRARY})
  endif()
endif(NEMO_FETCH_DEPS)

if(spdlog_VERSION VERSION_GREATER_EQUAL "1.6.0")
  set(NEMO_SPDLOG_FWD ON)
endif()

# Include optional subdirectories
include_directories(
    SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/source
)

# Generate nemo_git_version.h in the build directory
execute_process(
  COMMAND git log --pretty=format:%h -n 1
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE NEMO_GIT_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
if(NOT NEMO_GIT_VERSION)
  set(NEMO_GIT_VERSION "unknown")
endif()

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/nemo_git_version.h.in
  ${PROJECT_BINARY_DIR}/generated/nemo_git_version.h
)

add_subdirectory(source)

if(NOT NEMO_FETCH_DEPS)
  # Install the nemo library, its headers, and the CMake package files so that
  # other projects can use find_package(NEMO) and link to NEMO::nemo
  include(CMakePackageConfigHelpers)
  set(NEMO_INSTALL_CMAKEDIR ${NEMO_INSTALL_LIBDIR}/cmake/NEMO)

  install(
      TARGETS nemo
      EXPORT NEMOTargets
      ARCHIVE DESTINATION ${NEMO_INSTALL_LIBDIR}
  )

  install(
      DIRECTORY ${PROJECT_SOURCE_DIR}/source/
      DESTINATION ${NEMO_INSTALL_INCLUDEDIR}/nemo
      FILES_MATCHING PATTERN "*.h"
  )

  install(
      EXPORT NEMOTargets
      NAMESPACE NEMO::
      DESTINATION ${NEMO_INSTALL_CMAKEDIR}
  )

  configure_package_config_file(
      ${PROJECT_SOURCE_DIR}/cmake/NEMOConfig.cmake.in
      ${PROJECT_BINARY_DIR}/NEMOConfig.cmake
      INSTALL_DESTINATION ${NEMO_INSTALL_CMAKEDIR}
  )

  write_basic_package_version_file(
      ${PROJECT_BINARY_DIR}/NEMOConfigVersion.cmake
      COMPATIBILITY SameMajorVersion
  )

  install(
      FILES
          ${PROJECT_BINARY_DIR}/NEMOConfig.cmake
          ${PROJECT_BINARY_DIR}/NEMOConfigVersion.cmake
      DESTINATION ${NEMO_INSTALL_CMAKEDIR}
  )
endif(NOT NEMO_FETCH_DEPS)

# Add unit tests to run:
if(NEMO_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(NEMO_PYTHON)
  add_subdirectory(python)
endif()

message(STATUS "Build options:")
message(STATUS "    CMAKE_BUILD_TYPE         = ${CMAKE_BUILD_TYPE}")
message(STATUS "    NEMO_OMP                 = ${NEMO_OMP}")
message(STATUS "    NEMO_FETCH_DEPS          = ${NEMO_FETCH_DEPS}")
message(STATUS "    NEMO_PYTHON              = ${NEMO_PYTHON}")
message(STATUS "    SUNDIALS version         = ${SUNDIALS_VERSION_STRING}")

