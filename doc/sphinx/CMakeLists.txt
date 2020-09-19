find_package(Sphinx REQUIRED)
find_package(Python COMPONENTS Interpreter Development)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(SPHINX_BUILD ${CMAKE_CURRENT_SOURCE_DIR}/build)

add_custom_target(Sphinx ALL
                  COMMAND
                  ${SPHINX_EXECUTABLE} -b html
		  -Dbreathe_projects.libecpint=${DOXYGEN_OUTPUT_DIR}
		  ${SPHINX_SOURCE} ${SPHINX_BUILD}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  COMMENT "Generating documentation with Sphinx")
set_target_properties(Sphinx PROPERTIES EXCLUDE_FROM_ALL TRUE)
