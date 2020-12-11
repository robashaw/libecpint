add_executable(type1_test2 test2.cpp)
target_link_libraries(type1_test2 ecpint)

add_test(NAME Type1Test2
	 COMMAND type1_test2
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/type1_test2)
