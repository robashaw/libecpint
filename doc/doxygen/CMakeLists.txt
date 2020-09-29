find_package(Doxygen)

# Find all the public headers
get_target_property(LIBECPINT_PUBLIC_HEADER_DIR ecpint  INTERFACE_INCLUDE_DIRECTORIES)
file(GLOB_RECURSE LIBECPINT_PUBLIC_HEADERS ${LIBECPINT_PUBLIC_HEADER_DIR}/*.hpp)

#This will be the main output of our command
set(DOXYGEN_INDEX_FILE ${CMAKE_CURRENT_SOURCE_DIR}/html/index.html)

add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS ${LIBECPINT_PUBLIC_HEADERS}
                   COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                   MAIN_DEPENDENCY Doxyfile
                   COMMENT "Generating docs")

add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})
set_target_properties(Doxygen PROPERTIES EXCLUDE_FROM_ALL TRUE)
