target_sources(nemo PRIVATE

  interpolate.cpp
  cvode_integrator.cpp
)

# Boost and SUNDIALS headers are included by NEMO's public headers
if(NOT NEMO_FETCH_DEPS)
  get_filename_component(SUNDIALS_TOP_INCLUDE_DIR ${SUNDIALS_INCLUDE_DIR} DIRECTORY)
  target_include_directories(nemo PUBLIC ${SUNDIALS_TOP_INCLUDE_DIR})
endif()
target_compile_definitions(nemo PUBLIC CVODE${SUNDIALS_VERSION_MAJOR})

target_link_libraries(nemo PUBLIC ${NEMO_BOOST_TARGET} ${NEMO_MPI_TARGET})
if(NEMO_FETCH_DEPS)
  # The fetched SUNDIALS headers come with the targets; PUBLIC so that
  # code including nemo.h finds them.
  target_link_libraries(nemo PUBLIC ${NEMO_SUNDIALS_TARGETS})
else()
  target_link_libraries(nemo PRIVATE ${NEMO_SUNDIALS_TARGETS})
endif()
