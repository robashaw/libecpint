# minimum version of CMake is 3.12
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# define project
project( libecpint
         VERSION 1.0.5
         LANGUAGES C CXX)

set(API_VERSION 1)
set (CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#build type
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type specified. Will build Release")
  set(CMAKE_BUILD_TYPE "RELEASE")
endif(NOT CMAKE_BUILD_TYPE)


# code generator variables
set (LIBECPINT_MAX_L "5" CACHE STRING "Maximum angular momentum")
set (LIBECPINT_MAX_UNROL "1" CACHE STRING "Maximum L for unrolling")

# configure the config header to pass the above variables to the program
configure_file (
  "${PROJECT_SOURCE_DIR}/include/libecpint/config.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/libecpint/config.hpp"
)

include(GNUInstallDirs)
include(ExternalProject)
include(external/ImportPugiXML.cmake)
message(${PUGIXML_INCLUDE_DIR})

add_subdirectory(external)
add_subdirectory(src)

option(LIBECPINT_BUILD_TESTS "Enables Libecpint tests." ON)
if (LIBECPINT_BUILD_TESTS)
  include(CTest)
  include(external/ImportGTest.cmake)
  enable_testing()
  add_subdirectory(tests)
endif()
option(LIBECPINT_BUILD_DOCS "Enables Libecpint documentation." ON)
if (LIBECPINT_BUILD_DOCS)
  add_subdirectory(doc)
endif()

install (DIRECTORY include/libecpint DESTINATION include)
install (FILES include/libecpint.hpp DESTINATION include)
install (FILES "${CMAKE_CURRENT_BINARY_DIR}/include/libecpint/config.hpp" DESTINATION include/libecpint)
install (DIRECTORY share/libecpint DESTINATION share)
