add_executable(bfield_test1 test1.cpp)
target_link_libraries(bfield_test1 ecpint)

add_test(NAME BFieldTest1
	 COMMAND bfield_test1
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/bfield_test1)
