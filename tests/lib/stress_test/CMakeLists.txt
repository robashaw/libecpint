add_executable(stress_test test.cpp)
target_link_libraries(stress_test ecpint)

add_custom_target(StressTest
		  DEPENDS stress_test)
		  
add_custom_command(TARGET StressTest
	           POST_BUILD
		   COMMAND stress_test ${CMAKE_SOURCE_DIR}/share/libecpint 
		  )
