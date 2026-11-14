find_package(Threads REQUIRED)
add_executable(thread_safety_test test.cpp)
target_link_libraries(thread_safety_test ecpint Threads::Threads)

add_test(NAME ThreadSafety
	 COMMAND thread_safety_test
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/thread_safety)
ecpint_test_requires_max_l(ThreadSafety 3)
