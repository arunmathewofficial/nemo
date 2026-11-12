add_executable(isochoric-nei isochoric-nei.cpp)
target_link_libraries(isochoric-nei PRIVATE nemo MPI::MPI_CXX)
target_compile_features(isochoric-nei PRIVATE cxx_std_14)
target_compile_definitions(isochoric-nei PRIVATE
    ISOCHORIC_NEI_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/../data/isochoric_nei_reference.txt"
)
if(NEMO_OMP)
    target_link_libraries(isochoric-nei PRIVATE OpenMP::OpenMP_CXX)
endif()

add_custom_target(test-isochoric-nei
    COMMAND $<TARGET_FILE:isochoric-nei>
    DEPENDS isochoric-nei
    VERBATIM
)
add_test(NAME isochoric-nei COMMAND isochoric-nei)
