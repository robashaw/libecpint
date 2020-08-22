add_executable(deriv_test2 test2.cpp)
target_link_libraries(deriv_test2 ecpint)

add_test(NAME DerivTest2
	 COMMAND deriv_test2
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/deriv_test2)
