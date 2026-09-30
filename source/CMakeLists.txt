
add_subdirectory(nemo)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(
  nemo STATIC

    nemo.cpp
    mellema_cooling.cpp
    atomic_physics_database.cpp
    coll_ionise_recomb.cpp
    chianti_cooling.cpp
    stellar_atmosphere_models.cpp
    multi_ion_photo_xsections.cpp
    photoionisation.cpp
    charge_exchange.cpp
)

target_link_libraries(global PRIVATE spdlog::spdlog)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    include(../cmake/development_flags.cmake)
    add_definitions("-DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE")
    target_link_libraries(global INTERFACE development_flags)
else()
    add_definitions("-DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_INFO")
endif(CMAKE_BUILD_TYPE MATCHES Debug)


add_executable(test-nemo test-nemo.cpp)
target_link_libraries(test-nemo PRIVATE nemo)
target_compile_features(test-nemo PRIVATE cxx_std_14)

install( 
    TARGETS test-nemo nemo
    EXPORT test-nemo-targets
    RUNTIME DESTINATION ${NEMO_INSTALL_BINDIR}
    LIBRARY DESTINATION ${NEMO_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${NEMO_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${NEMO_INSTALL_INCLUDEDIR}
)


