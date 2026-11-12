add_executable(radial_small_argument_test test.cpp)
target_link_libraries(radial_small_argument_test ecpint)

add_test(NAME RadialSmallArgument
	 COMMAND radial_small_argument_test
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/radial_small_argument)
