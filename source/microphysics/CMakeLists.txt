add_library(tools
  
  interpolate.cpp
)

target_link_libraries(tools PRIVATE Boost::boost)
target_link_libraries(tools PUBLIC spdlog::spdlog)

if (CMAKE_BUILD_TYPE MATCHES "Debug")
    message(STATUS "Curses_FOUND: ${Curses_FOUND}")
    if (Curses_FOUND)
        target_link_libraries(tools PRIVATE ${CURSES_LIBRARIES})
    endif(Curses_FOUND)
    target_link_libraries(tools PRIVATE ${READLINE_LIBRARY})
endif(CMAKE_BUILD_TYPE MATCHES "Debug")

add_library(
  cvode_int

  cvode_integrator.cpp
)

target_link_libraries(cvode_int PRIVATE tools ${SUNDIALS_LIBRARY})
if (${SUNDIALS_VERSION_STRING} VERSION_GREATER_EQUAL 7.0.0)
  target_link_libraries(cvode_int PRIVATE tools ${SUNDIALS_sundials_core_LIBRARY})
endif()
target_link_libraries(cvode_int PRIVATE spdlog::spdlog)


