include(FindPackageHandleStandardArgs)

if(NOT units_FOUND)
    find_path( units_INCLUDE_DIR
        HINTS /
            /usr/
            /Users/jm/.local/
            ${CMAKE_INSTALL_PREFIX}/
            ${CMAKE_INSTALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${units_DIR}/
        PATH_SUFFIXES include
                    include/units
                    include/x86_64-linux-gnu
                    include/x86_64-linux-gnu/units
        NAMES units.hpp
        DOC "units headers"
    )
    find_library( units_LIBRARY
        HINTS /
            /usr/
            /Users/jm/.local/
            ${CMAKE_INSTALL_PREFIX}/
            ${CMAKE_INSTALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${units_DIR}/
        PATH_SUFFIXES lib lib64
                    lib/units lib64/units
                    lib/x86_64-linux-gnu lib64/x86_64-linux-gnu
                    lib/x86_64-linux-gnu/units lib64/x86_64-linux-gnu/units
        NAMES libunits.a
        DOC "units libraries"
    )

    find_package_handle_standard_args(units
            FOUND_VAR units_FOUND
            REQUIRED_VARS units_INCLUDE_DIR units_LIBRARY units_VERSION_MAJOR
            VERSION_VAR units_VERSION_STRING
    )            
endif(NOT units_FOUND)
