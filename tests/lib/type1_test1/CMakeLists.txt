add_executable(type1_test1 test1.cpp)
target_link_libraries(type1_test1 ecpint)

add_test(NAME Type1Test1
	 COMMAND type1_test1
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/type1_test1)
