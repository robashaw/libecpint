add_executable(quad_test
	gaussquad.cpp
)

target_link_libraries(quad_test gtest)
target_link_libraries(quad_test ecpint)

add_test(
	NAME
	  GaussQuad
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/gaussquad/quad_test
)
