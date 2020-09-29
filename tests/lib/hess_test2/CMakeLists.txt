add_executable(hess_test2 test2.cpp)
target_link_libraries(hess_test2 ecpint)

add_test(NAME HessTest2
	 COMMAND hess_test2
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/hess_test2)
