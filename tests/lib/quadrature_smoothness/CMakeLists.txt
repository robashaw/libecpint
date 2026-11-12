add_executable(quadrature_smoothness_test test.cpp)
target_link_libraries(quadrature_smoothness_test ecpint)

add_test(NAME QuadratureSmoothness
	 COMMAND quadrature_smoothness_test
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/quadrature_smoothness)
