add_executable(cie-test cie-test.cpp)
target_link_libraries(cie-test PRIVATE nemo)
target_link_libraries(cie-test PRIVATE MPI::MPI_CXX)
target_compile_features(cie-test PRIVATE cxx_std_14)
if(NEMO_OMP)
    target_link_libraries(cie-test PRIVATE OpenMP::OpenMP_CXX)
endif(NEMO_OMP)

add_custom_target(
  test-cie-test
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run-cie-test.sh"
  VERBATIM
)
add_test(NAME cie-test COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run-cie-test.sh")


