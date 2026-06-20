add_executable(screen_scan_test test.cpp)
target_link_libraries(screen_scan_test ecpint)

add_test(NAME ScreenScan
	 COMMAND screen_scan_test
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/screen_scan)
