include(FindPackageHandleStandardArgs)

function(find_sundials_version)
    include(CheckCXXSymbolExists)
    include(CMakePushCheckState)
    set(CMAKE_REQUIRED_QUIET ON)

    CMAKE_PUSH_CHECK_STATE()

    set(CMAKE_REQUIRED_INCLUDES "${SUNDIALS_INCLUDE_DIR}")
    file(READ "${SUNDIALS_INCLUDE_DIR}/sundials_config.h" version_header)
    string(REGEX MATCH "#define SUNDIALS_VERSION \"[0-9]+.[0-9]+.[0-9]+\"" version_string_macrodef "${version_header}")
    string(REGEX MATCH "[0-9]+.[0-9]+.[0-9]+" SUNDIALS_VERSION_STRING "${version_string_macrodef}")
    string(REGEX MATCH "#define SUNDIALS_VERSION_MAJOR [0-9]+" version_major_macrodef "${version_header}")
    string(REGEX MATCH "[0-9]+" SUNDIALS_VERSION_MAJOR "${version_major_macrodef}")

    set(SUNDIALS_VERSION_MAJOR ${SUNDIALS_VERSION_MAJOR} PARENT_SCOPE)
    set(SUNDIALS_VERSION_STRING ${SUNDIALS_VERSION_STRING} PARENT_SCOPE)

    CMAKE_POP_CHECK_STATE()
endfunction()

if(NOT SUNDIALS_FOUND)
    find_path( SUNDIALS_INCLUDE_DIR
        HINTS /
            /usr/
            ${CMAKE_INSTALL_PREFIX}/
            ${CMAKE_INSTALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${SUNDIALS_DIR}/
        PATH_SUFFIXES include
                    include/sundials
                    include/x86_64-linux-gnu
                    include/x86_64-linux-gnu/sundials
        NAMES sundials_config.h
        DOC "Sundials headers"
    )
    find_sundials_version()    
    find_library( SUNDIALS_LIBRARY
        HINTS /
            /usr/
            ${CMAKE_INSTALL_PREFIX}/
            ${CMAKE_INSTALL_PREFIX}/../
            ${CMAKE_PREFIX_PATH}/
            ${SUNDIALS_DIR}/
        PATH_SUFFIXES lib lib64
                    lib/sundials lib64/sundials
                    lib/x86_64-linux-gnu lib64/x86_64-linux-gnu
                    lib/x86_64-linux-gnu/sundials lib64/x86_64-linux-gnu/sundials
        NAMES libsundials_cvode.a libsundials_cvode.so
        DOC "Sundials libraries"
    )

    
# find the SUNDIALS libraries (from https://github.com/cadet/)
    set(SUNDIALS_WANT_COMPONENTS
    sundials_core
    sundials_cvode
    sundials_nvecserial
    )
    foreach(LIB ${SUNDIALS_WANT_COMPONENTS})
        if (UNIX AND SUNDIALS_PREFER_STATIC_LIBRARIES)
            # According to bug 1643 on the CMake bug tracker, this is the
            # preferred method for searching for a static library.
            # See http://www.cmake.org/Bug/view.php?id=1643.  We search
            # first for the full static library name, but fall back to a
            # generic search on the name if the static search fails.
            set(THIS_LIBRARY_SEARCH lib${LIB}.a ${LIB})
        else()
            set(THIS_LIBRARY_SEARCH ${LIB})
        endif()

        find_library(SUNDIALS_${LIB}_LIBRARY
            NAMES ${THIS_LIBRARY_SEARCH}
            PATHS
                ${SUNDIALS_USER_PATHS}
                ${SUNDIALS_ROOT}
                ${SUNDIALS_ROOT_DIR}
                ${PKGCONFIG_SUNDIALS_INCLUDE_DIRS}
                ${SUNDIALS_INCLUDE_DIR}/..
            ENV
                SUNDIALS_ROOT
                SUNDIALS_ROOT_DIR
            PATH_SUFFIXES
                lib
                Lib
        )

        set(SUNDIALS_${LIB}_FOUND FALSE)
        if (SUNDIALS_${LIB}_LIBRARY)
            list(APPEND SUNDIALS_LIBRARIES ${SUNDIALS_${LIB}_LIBRARY})
            set(SUNDIALS_${LIB}_FOUND TRUE)
        endif()
        mark_as_advanced(SUNDIALS_${LIB}_LIBRARY)
    endforeach()

    foreach(LIB ${SUNDIALS_WANT_COMPONENTS})
        if (SUNDIALS_${LIB}_LIBRARY)
            add_library(SUNDIALS::${LIB} UNKNOWN IMPORTED)
            set_target_properties(SUNDIALS::${LIB} PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${SUNDIALS_INCLUDE_DIRS}"
                IMPORTED_LOCATION ${SUNDIALS_${LIB}_LIBRARY}
            )
        endif()
    endforeach()

    find_package_handle_standard_args(SUNDIALS
            FOUND_VAR SUNDIALS_FOUND
            REQUIRED_VARS SUNDIALS_INCLUDE_DIR SUNDIALS_LIBRARY SUNDIALS_VERSION_MAJOR
            VERSION_VAR SUNDIALS_VERSION_STRING
    )            
endif(NOT SUNDIALS_FOUND)
