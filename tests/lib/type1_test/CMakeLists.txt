add_executable(type1_test test1.cpp)
target_link_libraries(type1_test ecpint)

add_test(NAME Type1Test
	 COMMAND type1_test
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/type1_test)
