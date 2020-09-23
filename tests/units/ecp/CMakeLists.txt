add_executable(ecp_test
	main.cpp
	ecp.cpp
	ecpbasis.cpp
)

message(${GTEST_INCLUDE_DIR})
target_link_libraries(ecp_test gtest)
target_link_libraries(ecp_test ecpint)

add_test(
	NAME
	  GaussianECP
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/ecp/ecp_test
	WORKING_DIRECTORY
	  ${CMAKE_SOURCE_DIR}/tests/units/ecp
)
