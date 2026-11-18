# `make docs` runs Doxygen (XML for breathe/exhale), then Sphinx. All output goes in the build tree,
# except exhale's generated rst, which it insists on putting in doc/sphinx/source/api (gitignored).
find_package(Doxygen)
find_package(Sphinx)

if (NOT DOXYGEN_FOUND OR NOT SPHINX_FOUND)
  set(_docs_msg "Cannot build docs:")
  if (NOT DOXYGEN_FOUND)
    string(APPEND _docs_msg " doxygen not found (set DOXYGEN_EXECUTABLE).")
  endif()
  if (NOT SPHINX_FOUND)
    string(APPEND _docs_msg " sphinx-build not found (set SPHINX_EXECUTABLE).")
  endif()
  string(APPEND _docs_msg " Python requirements are in ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/requirements.txt.")
  message(STATUS "${_docs_msg}")
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/docs_unavailable.cmake "message(FATAL_ERROR \"${_docs_msg}\")\n")
  add_custom_target(docs COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/docs_unavailable.cmake)
  return()
endif()

# Doxygen: wrap the checked-in Doxyfile so that paths are absolute and output lands here. Only the
# installed headers are documented: reading the sources too makes Doxygen merge declarations with
# their definitions, whose signatures differ in top-level const, and breathe then renders nothing.
set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUT}/xml/index.xml)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  "@INCLUDE = ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile\n"
  "OUTPUT_DIRECTORY = ${DOXYGEN_OUT}\n"
  "INPUT = ${PROJECT_SOURCE_DIR}/README.md ${PROJECT_SOURCE_DIR}/include/libecpint.hpp ${PROJECT_SOURCE_DIR}/include/libecpint\n"
  "USE_MDFILE_AS_MAINPAGE = ${PROJECT_SOURCE_DIR}/README.md\n"
  "STRIP_FROM_PATH = ${PROJECT_SOURCE_DIR}\n")

file(GLOB_RECURSE LIBECPINT_DOC_SOURCES CONFIGURE_DEPENDS
     ${PROJECT_SOURCE_DIR}/include/libecpint.hpp ${PROJECT_SOURCE_DIR}/include/libecpint/*.hpp)

add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile ${PROJECT_SOURCE_DIR}/README.md
                           ${LIBECPINT_DOC_SOURCES}
                   COMMENT "Generating API docs with Doxygen")

# Sphinx: always rerun (it is incremental itself)
set(SPHINX_OUT ${CMAKE_CURRENT_BINARY_DIR}/html)
add_custom_target(docs
                  COMMAND ${SPHINX_EXECUTABLE} -b html
                          -Dbreathe_projects.libecpint=${DOXYGEN_OUT}/xml
                          ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/source ${SPHINX_OUT}
                  COMMAND ${CMAKE_COMMAND} -E echo "Documentation written to ${SPHINX_OUT}/index.html"
                  DEPENDS ${DOXYGEN_INDEX_FILE}
                  COMMENT "Generating documentation with Sphinx")
