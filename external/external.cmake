add_custom_target(External)

add_subdirectory(external/Faddeeva)
include_directories(${CMAKE_SOURCE_DIR}/external/Faddeeva)
set(FADDEEVA_LIBRARY ${CMAKE_BINARY_DIR}/external/Faddeeva/libFaddeeva${CMAKE_STATIC_LIBRARY_SUFFIX})

set(EXTERNAL_SOURCE_DIR ${PROJECT_SOURCE_DIR}/external)
set(EXTERNAL_BUILD_DIR ${PROJECT_BINARY_DIR}/external/build)

#
# add pugixml

find_library(PUGIXML_LIBRARY pugixml
					HINTS "/usr/local" "/usr/local/lib" "/usr/lib")
find_path(PUGIXML_INCLUDE_DIR pugixml.hpp
					HINTS "/usr/local/include" "/usr/include")

if((NOT PUGIXML_LIBRARY) OR (NOT PUGIXML_INCLUDE_DIR))
	message("Unable to find pugixml, cloning...")
	
	ExternalProject_Add(pugixml_external
			PREFIX ${EXTERNAL_BUILD_DIR}/pugixml
			GIT_REPOSITORY https://github.com/zeux/pugixml
			CMAKE_ARGS “-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>”
			UPDATE_COMMAND ""
			INSTALL_COMMAND ""
			)
									       
	set(PUGIXML_LIBRARY ${EXTERNAL_BUILD_DIR}/pugixml/src/pugixml_external-build/libpugixml${CMAKE_STATIC_LIBRARY_SUFFIX})
	set(PUGIXML_INCLUDE_DIR ${EXTERNAL_BUILD_DIR}/pugixml/src/pugixml_external/src)
	
	add_library(libpugixml STATIC IMPORTED)
	set_target_properties(libpugixml PROPERTIES IMPORTED_LOCATION ${PUGIXML_LIBRARY})
									       
else()
	message("Found pugixml")
	
	add_library(libpugixml INTERFACE)
	target_include_directories(libpugixml INTERFACE ${PUGIXML_INCLUDE_DIR})
	target_link_libraries(libpugixml INTERFACE ${PUGIXML_LIBRARY})
endif()


#
# add google test

find_library(GTEST_LIBRARY gtest
					HINTS "/usr/local" "/usr/local/lib" "/usr/lib")
find_path(GTEST_INCLUDE_DIR gtest/gtest.h
					HINTS "/usr/local/include" "/usr/include")

if((NOT GTEST_LIBRARY) OR (NOT GTEST_INCLUDE_DIR))
	message("Unable to find google test, cloning...")
	
	# Download and unpack googletest at configure time
	configure_file(${CMAKE_SOURCE_DIR}/external/CMakeLists.txt.in googletest-download/CMakeLists.txt)
	execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
	  RESULT_VARIABLE result
	  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
	if(result)
	  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
	endif()
	execute_process(COMMAND ${CMAKE_COMMAND} --build .
	  RESULT_VARIABLE result
	  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
	if(result)
	  message(FATAL_ERROR "Build step for googletest failed: ${result}")
	endif()

	# Prevent overriding the parent project's compiler/linker
	# settings on Windows
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

	# Add googletest directly to our build. This defines
	# the gtest and gtest_main targets.
	add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
	                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
	                 )
					 
	set(GTEST_INCLUDE_DIR ${gtest_SOURCE_DIR}/include)
									       
else()
	message("Found googletest")
	
	add_library(gtest INTERFACE)
	target_include_directories(gtest INTERFACE ${GTEST_INCLUDE_DIR})
	target_link_libraries(gtest INTERFACE ${GTEST_LIBRARY})
endif()
