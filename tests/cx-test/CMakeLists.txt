add_executable(cx-test cx-test.cpp)
target_link_libraries(cx-test PRIVATE nemo tools cvode_int constants)
target_link_libraries(cx-test PRIVATE MPI::MPI_CXX)
target_compile_features(cx-test PRIVATE cxx_std_14)
if(NEMO_OMP)
    target_link_libraries(cx-test PRIVATE OpenMP::OpenMP_CXX)
endif(NEMO_OMP)

install(
    TARGETS cx-test
    EXPORT cx-test-targets
    RUNTIME DESTINATION ${NEMO_INSTALL_BINDIR}
    LIBRARY DESTINATION ${NEMO_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${NEMO_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${NEMO_INSTALL_INCLUDEDIR}
)

add_custom_target(
  test-cx-test
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run-cx-test.sh"
  VERBATIM
)
add_test(NAME cx-test COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run-cx-test.sh")


