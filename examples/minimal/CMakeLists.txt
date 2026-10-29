cmake_minimum_required(VERSION 3.15)
project(nemo_example LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

# Finds NEMO's installed package files (lib/cmake/NEMO/NEMOConfig.cmake).
# Point CMake at the install prefix with e.g. -DCMAKE_PREFIX_PATH=$HOME/.local
find_package(NEMO 1.0 REQUIRED)

add_executable(nemo-example main.cpp)

# Brings in libnemo.a plus all of NEMO's include paths, compile definitions
# and dependencies (SUNDIALS, Boost, MPI, spdlog).
target_link_libraries(nemo-example PRIVATE NEMO::nemo)
