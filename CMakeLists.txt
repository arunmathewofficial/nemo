cmake_minimum_required(VERSION 3.15.0)

project(
    NEMO
    VERSION 1.0.0
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


#  This is required on OSX so that the Boost code will compile successfully
set(CMAKE_CXX_STANDARD 14)

# Find dependant libraries
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

find_package(Boost 1.73.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(SUNDIALS REQUIRED)
include_directories(${SUNDIALS_INCLUDE_DIR}/../)
add_definitions("-DCVODE${SUNDIALS_VERSION_MAJOR}")
#message(STATUS "SUNDIALS version = ${SUNDIALS_VERSION_MAJOR}.${SUNDIALS_VERSION_MINOR}")
#message(STATUS "SUNDIALS include = ${SUNDIALS_INCLUDE_DIR}/../")

find_package(MPI REQUIRED)  # Homebrew sundials needs open MPI

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    find_package(Curses)
    find_package(READLINE REQUIRED)
    include_directories(${READLINE_INCLUDE_DIR})
endif(CMAKE_BUILD_TYPE MATCHES "Debug")

find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)

add_definitions("-DSPDLOG_FMT_EXTERNAL")
add_definitions("-DSPDLOG_COMPILED_LIB")

if (${spdlog_VERSION_MAJOR} GREATER_EQUAL 1)
    if (${spdlog_VERSION_MINOR} GREATER_EQUAL 6)
        message(STATUS "spdlog version = ${spdlog_VERSION_MAJOR}.${spdlog_VERSION_MINOR}")
        add_definitions("-DSPDLOG_FWD")
    endif (${spdlog_VERSION_MINOR} GREATER_EQUAL 6)
endif (${spdlog_VERSION_MAJOR} GREATER_EQUAL 1)

# Include optional subdirectories
include_directories(
    SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/source
)

add_subdirectory(source)
add_subdirectory(tests)

execute_process(
  COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tools/get-git-ver.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
)


message(STATUS "Build options:")
message(STATUS "    CMAKE_BUILD_TYPE         = ${CMAKE_BUILD_TYPE}")
