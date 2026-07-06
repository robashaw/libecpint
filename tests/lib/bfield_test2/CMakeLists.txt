add_executable(bfield_test2 test2.cpp)
target_link_libraries(bfield_test2 ecpint)

add_test(NAME BFieldTest2
	 COMMAND bfield_test2
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/bfield_test2)
