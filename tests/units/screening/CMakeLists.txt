add_executable(screening_test
	screening.cpp
)

target_link_libraries(screening_test ecpint_testsupport)

add_test(
	NAME
	  Screening
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/screening/screening_test
)
