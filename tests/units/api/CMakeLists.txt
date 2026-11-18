add_executable(api_test
	api.cpp
)

target_link_libraries(api_test ecpint_testsupport)

add_test(
	NAME
	  Api
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/api/api_test
)
