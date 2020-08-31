# minimum version of CMake is 3.2
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

# define project
project( libecpint
         VERSION 0.9
         LANGUAGES C CXX)

# we need at least C++11
set (CMAKE_CXX_STANDARD 11)

# code generator variables
set (LIBECPINT_MAX_L "5" CACHE STRING "Maximum angular momentum")
set (LIBECPINT_MAX_UNROL "2" CACHE STRING "Maximum L for unrolling")

# configure the config header to pass the above variables to the program
configure_file (
	      "${PROJECT_SOURCE_DIR}/include/libecpint/config.hpp.in"
	      "${CMAKE_CURRENT_BINARY_DIR}/include/libecpint/config.hpp"
)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/libecpint)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/libecpint)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/Faddeeva)

add_subdirectory(src)

enable_testing()
add_subdirectory(tests)

install (DIRECTORY include/libecpint DESTINATION include)
install (FILES include/libecpint.hpp DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/include/libecpint/config.hpp" DESTINATION include/libecpint)
