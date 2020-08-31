execute_process(COMMAND python "${CMAKE_CURRENT_SOURCE_DIR}/makelist.py" "${LIBECPINT_MAX_L}" "${CMAKE_CURRENT_BINARY_DIR}")

file(READ "${CMAKE_CURRENT_BINARY_DIR}/qlist.txt" GENERATED_SOURCES)
STRING(REGEX REPLACE "\n" ";" GENERATED_SOURCES "${GENERATED_SOURCES}") 

add_executable(generate generate.cpp 
			lib/mathutil.cpp
			lib/angular.cpp)

execute_process(COMMAND mkdir -p "${CMAKE_CURRENT_BINARY_DIR}/generated")
execute_process(COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/generated/ecpint_gen.part" "${CMAKE_CURRENT_BINARY_DIR}/generated")
execute_process(COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/generated/qgen.part" "${CMAKE_CURRENT_BINARY_DIR}/generated")
add_custom_command(
	OUTPUT ${GENERATED_SOURCES}
	COMMAND generate "${PROJECT_SOURCE_DIR}/include/libecpint/"
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/generated/ecpint_gen.part 
	COMMENT "Generating and unrolling integral code."
)

file(GLOB_RECURSE LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cpp")

add_library(ecpint ${LIB_SOURCES} ${GENERATED_SOURCES})
if(NOT BUILD_SHARED_LIBS)
  set_property(TARGET ecpint PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
target_include_directories(ecpint
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/libecpint>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/libecpint>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
) 

install (TARGETS ecpint DESTINATION lib)

