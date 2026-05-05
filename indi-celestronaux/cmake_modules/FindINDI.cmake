find_path(INDI_INCLUDE_DIR NAMES indiapi.h PATH_SUFFIXES libindi)
find_library(INDI_LIBRARIES NAMES indidriver)
find_library(INDI_ALIGNMENT_LIBRARIES NAMES indiAlignmentDriver)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(INDI DEFAULT_MSG INDI_LIBRARIES INDI_INCLUDE_DIR)

if(INDI_FOUND)
    set(INDI_INCLUDE_DIRS ${INDI_INCLUDE_DIR})
    set(INDI_LIBRARIES ${INDI_LIBRARIES} ${INDI_ALIGNMENT_LIBRARIES})
    if(NOT TARGET INDI::INDI)
        add_library(INDI::INDI UNKNOWN IMPORTED)
        set_target_properties(INDI::INDI PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${INDI_INCLUDE_DIRS}"
            IMPORTED_LOCATION "${INDI_LIBRARIES}"
        )
    endif()
endif()
