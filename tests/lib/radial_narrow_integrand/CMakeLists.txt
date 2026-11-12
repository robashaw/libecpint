add_executable(radial_narrow_integrand_test test.cpp)
target_link_libraries(radial_narrow_integrand_test ecpint)

add_test(NAME RadialNarrowIntegrand
	 COMMAND radial_narrow_integrand_test
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/radial_narrow_integrand)
