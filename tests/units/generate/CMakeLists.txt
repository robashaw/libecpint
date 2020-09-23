add_executable(gen_test
	generate.cpp
)

target_link_libraries(gen_test gtest)
target_link_libraries(gen_test ecpint)

add_test(
	NAME
	  Generator
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/generate/gen_test
)
