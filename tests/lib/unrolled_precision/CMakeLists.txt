add_executable(unrolled_precision_test test.cpp)
target_link_libraries(unrolled_precision_test ecpint)

add_test(NAME UnrolledPrecision
	 COMMAND unrolled_precision_test
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/unrolled_precision)
