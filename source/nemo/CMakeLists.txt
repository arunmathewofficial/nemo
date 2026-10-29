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
add_library(NEMO::nemo ALIAS nemo)

target_link_libraries(nemo PUBLIC spdlog::spdlog)
if(NEMO_SPDLOG_FWD)
  target_compile_definitions(nemo PUBLIC SPDLOG_FWD)
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
    target_link_libraries(nemo INTERFACE $<BUILD_INTERFACE:development_flags>)
endif(CMAKE_BUILD_TYPE MATCHES Debug)


target_include_directories(nemo PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>
  $<INSTALL_INTERFACE:${NEMO_INSTALL_INCLUDEDIR}/nemo>
)
# Generated nemo_git_version.h (not installed)
target_include_directories(nemo PRIVATE ${PROJECT_BINARY_DIR}/generated)

