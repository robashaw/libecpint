add_executable(test2 test2.cpp)
target_link_libraries(test2 ecpint)

add_test(NAME IntTest2
	 COMMAND test2 
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/int_test2)
