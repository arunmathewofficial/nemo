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

