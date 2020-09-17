# minimum version of CMake is 3.2
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

# define project
project( libecpint
         VERSION 0.9
         LANGUAGES C CXX)

# code generator variables
set (LIBECPINT_MAX_L "5" CACHE STRING  "Maximum angular momentum")
set (LIBECPINT_MAX_UNROL "1" CACHE STRING "Maximum L for unrolling")

# configure the config header to pass the above variables to the program
configure_file (
	      "${PROJECT_SOURCE_DIR}/include/libecpint/config.hpp.in"
	      "${CMAKE_CURRENT_BINARY_DIR}/include/libecpint/config.hpp"
)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/libecpint)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/libecpint)

add_subdirectory(src)

option(LIBECPINT_BUILD_TESTS "Enables Libecpint tests." ON)
if (LIBECPINT_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

install(DIRECTORY include/libecpint DESTINATION include)
install(FILES include/libecpint.hpp DESTINATION include)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/libecpint/config.hpp DESTINATION include/libecpint)
