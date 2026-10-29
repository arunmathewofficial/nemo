
add_subdirectory(nemo)
add_subdirectory(solver)

target_sources(nemo PRIVATE

  nemo-constants.cpp
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(CMAKE_BUILD_TYPE MATCHES Debug)
    include(../cmake/development_flags.cmake)
    add_definitions("-DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE")
    target_link_libraries(nemo INTERFACE $<BUILD_INTERFACE:development_flags>)
else()
    add_definitions("-DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_INFO")
endif(CMAKE_BUILD_TYPE MATCHES Debug)
