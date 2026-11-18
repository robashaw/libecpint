add_executable(assembly_test
	assembly.cpp
)

target_link_libraries(assembly_test gtest)
target_link_libraries(assembly_test ecpint)

add_test(
	NAME
	  Assembly
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/assembly/assembly_test
)
