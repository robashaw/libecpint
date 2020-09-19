set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/xml)

add_subdirectory(doxygen)
add_subdirectory(sphinx)

add_custom_target(docs DEPENDS Doxygen Sphinx)
set_target_properties(docs PROPERTIES EXCLUDE_FROM_ALL TRUE)
