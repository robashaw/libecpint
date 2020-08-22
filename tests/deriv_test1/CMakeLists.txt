add_executable(deriv_test1 test1.cpp)
target_link_libraries(deriv_test1 ecpint)

add_test(NAME DerivTest1
	 COMMAND deriv_test1
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/deriv_test1)
