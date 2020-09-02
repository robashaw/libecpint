# minimum version of CMake is 3.2
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

# define project
project( libecpint
         VERSION 0.9
         LANGUAGES C CXX)

# code generator variables
option (LIBECPINT_MAX_L "Maximum angular momentum" "5")
option (LIBECPINT_MAX_UNROL "Maximum L for unrolling" "2")

# configure the config header to pass the above variables to the program
configure_file (
	      "${PROJECT_SOURCE_DIR}/include/libecpint/config.hpp.in"
	      "${PROJECT_SOURCE_DIR}/include/libecpint/config.hpp"
)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/libecpint)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/Faddeeva)

add_subdirectory(src)

option(LIBECPINT_BUILD_TESTS "Enables Libecpint tests." ON)
if (LIBECPINT_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

install (DIRECTORY include/libecpint DESTINATION include)
install (FILES include/libecpint.hpp DESTINATION include)
