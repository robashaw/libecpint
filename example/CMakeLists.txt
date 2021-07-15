cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# define project
project( LibecpintExample
         VERSION 1.0
         LANGUAGES CXX)

add_executable(example example.cpp)
target_link_libraries(example PRIVATE ecpint Faddeeva pugixml)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/iodine.bas
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/hydrogen.bas
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

