function(gather_sources FILEPATH_DIRS_ARG CIMEROOT_ARG)
  set(SRCROOT_REL "${CIMEROOT_ARG}/..")
  set(BASENAME_SET)
  set(SOURCES_RESULT)
  set(GEN_F90_SOURCES_RESULT)

  file(TO_CMAKE_PATH ${SRCROOT_REL} SRCROOT_ABS)
  foreach(DIRSEARCH ${FILEPATH_DIRS_ARG})
    file(GLOB MATCHES RELATIVE "${SRCROOT_ABS}/components" "${DIRSEARCH}/*.[Ffc]" "${DIRSEARCH}/*.[Ff]90" "${DIRSEARCH}/*.cpp" "${DIRSEARCH}/*.F90.in")
    if (MATCHES)
      foreach (MATCH IN LISTS MATCHES)
        get_filename_component(BASENAME ${MATCH} NAME)
        list(FIND BASENAME_SET ${BASENAME} BASENAME_WAS_FOUND)
        if (BASENAME_WAS_FOUND EQUAL -1)
          list(APPEND SOURCES_RESULT ${MATCH})
          list(APPEND BASENAME_SET ${BASENAME})
        else()
          message(WARNING "Skipping repeated base filename ${BASENAME} for ${MATCH}")
        endif()
      endforeach()
    endif()
  endforeach()

  foreach(SOURCE_FILE IN LISTS SOURCES_RESULT)
    get_filename_component(SOURCE_EXT ${SOURCE_FILE} EXT)
    if (SOURCE_EXT STREQUAL ".F90.in")
      string(REPLACE ".in" "" SOURCE_NO_IN ${SOURCE_FILE})
      list(APPEND GEN_F90_SOURCES_RESULT ${SOURCE_NO_IN})
      list(APPEND SOURCES_RESULT ${SOURCE_NO_IN})
      list(REMOVE_ITEM SOURCES_RESULT ${SOURCE_FILE})
    endif()
  endforeach()

  # Return data to parent
  set(SOURCES_RESULT ${SOURCES_RESULT} PARENT_SCOPE)
  set(GEN_F90_SOURCES_RESULT ${GEN_F90_SOURCES_RESULT} PARENT_SCOPE)

endfunction()