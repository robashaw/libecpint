add_executable(analytic_test
	analytic.cpp
)

target_link_libraries(analytic_test gtest)
target_link_libraries(analytic_test ecpint_testsupport)
target_compile_definitions(analytic_test PRIVATE ANALYTIC_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_test(
	NAME
	  Analytic
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/analytic/analytic_test
)
