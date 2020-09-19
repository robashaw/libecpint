add_executable(api_test1 test1.cpp)
target_link_libraries(api_test1 ecpint)

add_test(NAME APITest1
	 COMMAND api_test1
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/api_test1)
