add_executable(test1 test1.cpp)
target_link_libraries(test1 ecpint)

add_test(NAME Test1
	 COMMAND test1
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/test1)
