add_executable(radial_fixed_quadrature_test test.cpp)
target_link_libraries(radial_fixed_quadrature_test ecpint)

add_test(NAME RadialFixedQuadrature
	 COMMAND radial_fixed_quadrature_test
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/radial_fixed_quadrature)
