add_executable(scaling_test test.cpp)
target_link_libraries(scaling_test ecpint)

add_custom_target(ScalingTest
		  DEPENDS scaling_test)
		  
add_custom_command(TARGET ScalingTest
	           POST_BUILD
		   COMMAND scaling_test ${CMAKE_SOURCE_DIR}/share/libecpint 
		  )