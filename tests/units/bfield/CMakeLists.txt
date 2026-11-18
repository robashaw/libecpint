add_executable(bfield_test
	bfield.cpp
)

target_link_libraries(bfield_test ecpint_testsupport)

add_test(
	NAME
	  BField
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/bfield/bfield_test
)
