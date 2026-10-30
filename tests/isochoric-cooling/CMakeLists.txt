add_executable(isochoric-cooling isochoric-cooling.cpp)
target_link_libraries(isochoric-cooling PRIVATE nemo MPI::MPI_CXX)
target_compile_features(isochoric-cooling PRIVATE cxx_std_14)
target_compile_definitions(isochoric-cooling PRIVATE
    ISOCHORIC_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/../data/isochoric_cooling_reference.txt"
)
if(NEMO_OMP)
    target_link_libraries(isochoric-cooling PRIVATE OpenMP::OpenMP_CXX)
endif()

add_custom_target(test-isochoric-cooling
    COMMAND $<TARGET_FILE:isochoric-cooling>
    DEPENDS isochoric-cooling
    VERBATIM
)
add_test(NAME isochoric-cooling COMMAND isochoric-cooling)
