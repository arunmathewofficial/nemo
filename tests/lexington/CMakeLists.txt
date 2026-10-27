
# HII40 Lexington test
add_executable(hii40-lexington hii40-lexington.cpp) # Build the executable
target_link_libraries(hii40-lexington PRIVATE
        nemo          # Microphysics solver
        tools         # Interpolation utilities
        cvode_int     # CVODE integration interface
        constants     # Physical constants
        common_functions  # some functions for testing code
        MPI::MPI_CXX  # MPI dependencies
)
target_compile_features(hii40-lexington PRIVATE cxx_std_14) # Require C++14
if(NEMO_OMP)
    target_link_libraries(hii40-lexington PRIVATE OpenMP::OpenMP_CXX)
endif(NEMO_OMP)

add_custom_target(
  test-hii40-lexington
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run-hii40-test.sh"
  VERBATIM
)
add_test(NAME hii40-lexington COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run-hii40-test.sh")


