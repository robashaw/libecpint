add_executable(gshell_test
	gshell.cpp
)

message(${GTEST_INCLUDE_DIR})
target_link_libraries(gshell_test gtest)
target_link_libraries(gshell_test ecpint)

add_test(
	NAME
	  GaussianShell
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/gshell/gshell_test
)
