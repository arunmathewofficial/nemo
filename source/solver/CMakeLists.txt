target_sources(nemo PRIVATE

  interpolate.cpp
  cvode_integrator.cpp
)

# Boost and SUNDIALS headers are included by NEMO's public headers
get_filename_component(SUNDIALS_TOP_INCLUDE_DIR ${SUNDIALS_INCLUDE_DIR} DIRECTORY)
target_include_directories(nemo PUBLIC ${SUNDIALS_TOP_INCLUDE_DIR})
target_compile_definitions(nemo PUBLIC CVODE${SUNDIALS_VERSION_MAJOR})

target_link_libraries(nemo PUBLIC Boost::boost MPI::MPI_CXX)
target_link_libraries(nemo PRIVATE ${SUNDIALS_LIBRARY})
if (${SUNDIALS_VERSION_STRING} VERSION_GREATER_EQUAL 7.0.0)
  target_link_libraries(nemo PRIVATE ${SUNDIALS_sundials_core_LIBRARY})
endif()
