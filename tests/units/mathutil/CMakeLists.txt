add_executable(math_test
	math.cpp
)

target_link_libraries(math_test gtest)
target_link_libraries(math_test ecpint)

add_test(
	NAME
	  MathUtil
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/mathutil/math_test
)
