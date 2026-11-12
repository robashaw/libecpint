add_executable(radial_closed_forms_test test.cpp)
target_link_libraries(radial_closed_forms_test ecpint)

add_test(NAME RadialClosedForms
	 COMMAND radial_closed_forms_test
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/radial_closed_forms)
