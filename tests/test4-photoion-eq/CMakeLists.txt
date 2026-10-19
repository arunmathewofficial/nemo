
# HII40 Lexington test
add_executable(photoionization_eq photoionization_eq.cpp) # Build the executable
target_link_libraries(photoionization_eq PRIVATE
        nemo          # Microphysics solver
        tools         # Interpolation utilities
        cvode_int     # CVODE integration interface
        constants     # Physical constants
        common_functions  # some functions for testing code
        MPI::MPI_CXX  # MPI dependencies
)
target_compile_features(photoionization_eq PRIVATE cxx_std_14) # Require C++14
if(NEMO_OMP)
    target_link_libraries(photoionization_eq PRIVATE OpenMP::OpenMP_CXX)
endif(NEMO_OMP)

add_custom_target(
  test-pieq
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run-pieq.sh"
  VERBATIM
)
add_test(NAME pieq COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run-pieq.sh")


