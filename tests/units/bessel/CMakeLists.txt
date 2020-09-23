add_executable(bessel_test
	bessel.cpp
)

target_link_libraries(bessel_test gtest)
target_link_libraries(bessel_test ecpint)

add_test(
	NAME
	  Bessel
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/bessel/bessel_test
)
