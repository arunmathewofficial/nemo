
add_subdirectory(nemo)
add_subdirectory(microphysics)

add_library(
  constants

  nemo-constants.cpp
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

target_link_libraries(nemo PRIVATE spdlog::spdlog)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    include(../cmake/development_flags.cmake)
    add_definitions("-DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE")
    target_link_libraries(nemo INTERFACE development_flags)
else()
    add_definitions("-DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_INFO")
endif(CMAKE_BUILD_TYPE MATCHES Debug)


add_executable(test-nemo test-nemo.cpp)
target_link_libraries(test-nemo PRIVATE nemo tools cvode_int constants)
target_link_libraries(test-nemo PRIVATE MPI::MPI_CXX)
target_compile_features(test-nemo PRIVATE cxx_std_14)

install( 
    TARGETS test-nemo nemo
    EXPORT test-nemo-targets
    RUNTIME DESTINATION ${NEMO_INSTALL_BINDIR}
    LIBRARY DESTINATION ${NEMO_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${NEMO_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${NEMO_INSTALL_INCLUDEDIR}
)


