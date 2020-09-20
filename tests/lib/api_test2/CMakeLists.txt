add_executable(api_test2 test2.cpp)
target_link_libraries(api_test2 ecpint)

add_test(NAME APITest2
	 COMMAND api_test2 ${CMAKE_SOURCE_DIR}/share/libecpint
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/api_test2)
