add_executable(support_test
	support.cpp
)

target_link_libraries(support_test ecpint_testsupport)

add_test(
	NAME
	  TestSupport
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/support/support_test
)
