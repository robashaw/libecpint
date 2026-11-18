cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# define project
project( LibecpintExample
         VERSION 1.0
         LANGUAGES CXX)

# Built standalone, use the installed libecpint; built with libecpint's tests, use its target
if (NOT TARGET ECPINT::ecpint)
  find_package(ecpint CONFIG REQUIRED)
endif()

add_executable(example example.cpp)
target_link_libraries(example PRIVATE ECPINT::ecpint)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/iodine.bas
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/hydrogen.bas
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME example
      COMMAND example ${CMAKE_SOURCE_DIR}/share/libecpint
	 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
