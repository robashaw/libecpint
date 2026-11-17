add_executable(derivatives_test
	derivatives.cpp
)

target_link_libraries(derivatives_test gtest)
target_link_libraries(derivatives_test ecpint)

add_test(
	NAME
	  DerivativeRecurrences
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/derivatives/derivatives_test
)
