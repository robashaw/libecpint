add_executable(radial_test
	radial.cpp
)

target_link_libraries(radial_test gtest)
target_link_libraries(radial_test ecpint)

# One ctest entry per test (Radial.Type2.LowOrders, ...), so that ctest -j runs them in parallel
include(GoogleTest)
gtest_discover_tests(radial_test TEST_PREFIX "Radial.")
