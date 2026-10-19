
# One-zone test
add_executable(one-zone one-zone.cpp)
target_link_libraries(one-zone PRIVATE nemo tools cvode_int constants)
target_link_libraries(one-zone PRIVATE MPI::MPI_CXX)
target_compile_features(one-zone PRIVATE cxx_std_14)
install( 
    TARGETS one-zone
    EXPORT one-zone-targets
    RUNTIME DESTINATION ${NEMO_INSTALL_BINDIR}
    LIBRARY DESTINATION ${NEMO_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${NEMO_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${NEMO_INSTALL_INCLUDEDIR}
)

# HII40 Lexington test
add_executable(hii40-lexington hii40-lexington.cpp) # Build the executable
target_link_libraries(hii40-lexington PRIVATE
        nemo          # Microphysics solver
        tools         # Interpolation utilities
        cvode_int     # CVODE integration interface
        constants     # Physical constants
        MPI::MPI_CXX  # MPI dependencies
)
target_compile_features(hii40-lexington PRIVATE cxx_std_14) # Require C++14
if(NEMO_OMP)
    target_link_libraries(hii40-lexington PRIVATE OpenMP::OpenMP_CXX)
endif(NEMO_OMP)

# CIE test
add_executable(cie-test cie-test.cpp)
target_link_libraries(cie-test PRIVATE nemo tools cvode_int constants)
target_link_libraries(cie-test PRIVATE MPI::MPI_CXX)
target_compile_features(cie-test PRIVATE cxx_std_14)
if(NEMO_OMP)
    target_link_libraries(cie-test PRIVATE OpenMP::OpenMP_CXX)
endif(NEMO_OMP)

install(
    TARGETS cie-test
    EXPORT cie-test-targets
    RUNTIME DESTINATION ${NEMO_INSTALL_BINDIR}
    LIBRARY DESTINATION ${NEMO_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${NEMO_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${NEMO_INSTALL_INCLUDEDIR}
)

# Add unit tests to run: 
enable_testing()

add_custom_target(
  test-one-zone
  COMMAND one-zone
  VERBATIM
)
add_test(NAME one-zone COMMAND one-zone)

add_custom_target(test-cie-test COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run-cie-test.sh" VERBATIM)
add_test(NAME cie-test COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run-cie-test.sh")

add_custom_target(test-hii40-lexington COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run-hii40-test.sh" VERBATIM)
add_test(NAME hii40-lexington COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run-hii40-test.sh")

# Glob all text files in the subdirectory
file(GLOB DATA_FILES "${CMAKE_CURRENT_SOURCE_DIR}/data/*.txt")
add_custom_target(copy_data ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data
        ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMENT "Copying data files..."
)
# Make your main target depend on it
add_dependencies(cie-test copy_data)
add_dependencies(hii40-lexington copy_data)
