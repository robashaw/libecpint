add_executable(shared_grid_test
	shared_grid.cpp
)

target_link_libraries(shared_grid_test gtest)
target_link_libraries(shared_grid_test ecpint)

add_test(
	NAME
	  SharedGrid
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/shared_grid/shared_grid_test
)
