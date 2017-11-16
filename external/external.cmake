add_custom_target(External)

add_subdirectory(external/Faddeeva)
include_directories(${CMAKE_SOURCE_DIR}/external/Faddeeva)
set(FADDEEVA_LIBRARY ${CMAKE_BINARY_DIR}/external/Faddeeva/libFaddeeva${CMAKE_STATIC_LIBRARY_SUFFIX})

set(EXTERNAL_SOURCE_DIR ${PROJECT_SOURCE_DIR}/external)
set(EXTERNAL_BUILD_DIR ${PROJECT_BINARY_DIR}/external/build)

#
# add pugixml
#
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
									       )
									       
									       set(PUGIXML_LIBRARY ${EXTERNAL_BUILD_DIR}/pugixml/src/pugixml_external-build/libpugixml${CMAKE_STATIC_LIBRARY_SUFFIX})
									       set(PUGIXML_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/pugixml/src/pugixml_external/include)
									       
									       add_library(libpugixml STATIC IMPORTED)
									       set_target_properties(libpugixml PROPERTIES IMPORTED_LOCATION ${PUGIXML_LIBRARY})
									       
else()
	message("Found pugixml")
	
	add_library(libpugixml INTERFACE)
	target_include_directories(libpugixml INTERFACE ${PUGIXML_INCLUDE_DIR})
	target_link_libraries(libpugixml INTERFACE ${PUGIXML_LIBRARY})
endif()
