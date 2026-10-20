
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
add_custom_target(
  test-one-zone
  COMMAND one-zone
  VERBATIM
)
add_test(NAME one-zone COMMAND one-zone)

# Add data files from sub-directory to cmake
file(GLOB DATA_FILES "${CMAKE_CURRENT_SOURCE_DIR}/data/*.txt")
add_custom_target(copy_data ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data
        ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMENT "Copying data files..."
)


add_library(
  common_functions

  common_functions.cpp
)



# CIE test
add_subdirectory(cie-test)
add_dependencies(cie-test copy_data)

# CX test
add_subdirectory(cx-test)
add_dependencies(cx-test copy_data)

# photoionization equilibrium test
add_subdirectory(test4-photoion-eq)
add_dependencies(photoionization_eq copy_data)

# HII40 Lexington test
add_subdirectory(lexington)
add_dependencies(hii40-lexington copy_data)