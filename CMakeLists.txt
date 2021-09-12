# minimum version of CMake is 3.12
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

include(CheckCXXCompilerFlag)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

#build type
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type specified. Will build Release")
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Release/Debug/RelWithDebInfo)")
else()
  string(TOUPPER "${CMAKE_BUILD_TYPE}" _upper_build_type)
  if("${_upper_build_type}" STREQUAL "DEBUG")
    add_compile_definitions(DEBUG)
  endif()
endif(NOT CMAKE_BUILD_TYPE)

# define project
project( libecpint
         VERSION 1.1.0
         LANGUAGES C CXX)

set(API_VERSION 1)
set (CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

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
option(LIBECPINT_USE_PUGIXML "Use pugixml and build file reading routines." ON)
if (LIBECPINT_USE_PUGIXML)
  include(external/ImportPugiXML.cmake)
  add_compile_definitions(HAS_PUGIXML)
else()
  message(STATUS "pugixml not configured; reading ECP definitions from files is disabled.")
endif()

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

find_package(CLANG_FORMAT 11)
if(CLANG_FORMAT_FOUND)
set(FORMAT_SOURCES)
list(APPEND FORMAT_SUBDIRS tests "include/libecpint" src)
foreach(DIR ${FORMAT_SUBDIRS})
    set(ABS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${DIR})
    file(GLOB_RECURSE ${DIR}_SOURCES ${ABS_DIR}/*.cpp ${ABS_DIR}/*.hpp)
    list(APPEND FORMAT_SOURCES ${${DIR}_SOURCES})
  endforeach()
    file(GLOB FORMAT_INCLUDES LIST_DIRECTORIES false ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
    list (APPEND FORMAT_SOURCES ${FORMAT_INCLUDES})
    add_custom_target(format
    COMMAND ${CLANG_FORMAT_EXECUTABLE} -i -style=file ${FORMAT_SOURCES}
    DEPENDS ${FORMAT_SOURCES})
endif()


install (DIRECTORY include/libecpint DESTINATION include)
install (FILES include/libecpint.hpp DESTINATION include)
install (FILES "${CMAKE_CURRENT_BINARY_DIR}/include/libecpint/config.hpp" DESTINATION include/libecpint)
install (DIRECTORY share/libecpint DESTINATION share)
