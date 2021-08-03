set(EXTERNAL_SOURCE_DIR ${PROJECT_SOURCE_DIR}/external)
set(EXTERNAL_BUILD_DIR ${PROJECT_BINARY_DIR}/external/build)

find_library(PUGIXML_LIBRARY pugixml
 HINTS "/usr/local" "/usr/local/lib" "/usr/lib"
)
find_path(PUGIXML_INCLUDE_DIR pugixml.hpp
 HINTS "/usr/local/include" "/usr/include"
)

if((NOT PUGIXML_LIBRARY) OR (NOT PUGIXML_INCLUDE_DIR))
  message(STATUS "Unable to find pugixml, cloning and building ...")

  ExternalProject_Add(pugixml_external
    PREFIX ${EXTERNAL_BUILD_DIR}/pugixml
    GIT_REPOSITORY https://github.com/zeux/pugixml
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>" "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    UPDATE_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${EXTERNAL_BUILD_DIR}/pugixml/src/pugixml_external-build/libpugixml${CMAKE_STATIC_LIBRARY_SUFFIX}
  )

  set(PUGIXML_LIBRARY ${EXTERNAL_BUILD_DIR}/pugixml/src/pugixml_external-build/libpugixml${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(PUGIXML_INCLUDE_DIR ${EXTERNAL_BUILD_DIR}/pugixml/src/pugixml_external/src)

  add_library(libpugixml STATIC IMPORTED)
  add_dependencies(libpugixml pugixml_external)
  set_target_properties(libpugixml PROPERTIES IMPORTED_LOCATION ${PUGIXML_LIBRARY})
  install(FILES ${PUGIXML_LIBRARY} DESTINATION lib)

else()
  message(STATUS "Found pugixml: ${PUGIXML_INCLUDE_DIR} ${PUGIXML_LIBRARY}")
  add_library(libpugixml INTERFACE)
  target_include_directories(libpugixml INTERFACE ${PUGIXML_INCLUDE_DIR})
  target_link_libraries(libpugixml INTERFACE ${PUGIXML_LIBRARY})
  install(TARGETS libpugixml EXPORT ecpintTargets DESTINATION lib)
endif()
