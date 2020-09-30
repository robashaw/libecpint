add_executable(multiarr_test
	multiarr.cpp
	twoindex.cpp
	threeindex.cpp
	fiveindex.cpp
	sevenindex.cpp
)

target_link_libraries(multiarr_test gtest)
target_link_libraries(multiarr_test ecpint)

add_test(
	NAME
	  MultiArray
	COMMAND
	  ${CMAKE_BINARY_DIR}/tests/units/multiarr/multiarr_test
)
