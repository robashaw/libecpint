add_executable(radial_high_order_test test.cpp)
target_link_libraries(radial_high_order_test ecpint)

add_test(NAME RadialHighOrder
	 COMMAND radial_high_order_test
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/radial_high_order)
