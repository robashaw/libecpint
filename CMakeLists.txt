# minimum version of CMake is 3.12
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

include(CheckCXXCompilerFlag)

# build type
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type specified. Will build Release.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Release/Debug/RelWithDebInfo)")
else()
  string(TOUPPER "${CMAKE_BUILD_TYPE}" _upper_build_type)
  if("${_upper_build_type}" STREQUAL "DEBUG")
    add_compile_definitions(DEBUG)
  endif()
endif(NOT CMAKE_BUILD_TYPE)

# define project
project( libecpint
         VERSION 1.1.0
         LANGUAGES CXX)

set(API_VERSION 1)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# code generator variables
set(LIBECPINT_MAX_UNROL "1" CACHE STRING "Maximum L for unrolling")

# per-thread cache for the radial factors of off-centre type 2 integrals, in MB
set(LIBECPINT_RADIAL_CACHE_MB "4" CACHE STRING "Per-thread radial factor cache size (MB)")
set(LIBECPINT_RADIAL_CACHE_MB_VALUES 1 2 4 8 16 32 64)
set_property(CACHE LIBECPINT_RADIAL_CACHE_MB PROPERTY STRINGS ${LIBECPINT_RADIAL_CACHE_MB_VALUES})
if (NOT LIBECPINT_RADIAL_CACHE_MB IN_LIST LIBECPINT_RADIAL_CACHE_MB_VALUES)
  message(FATAL_ERROR "LIBECPINT_RADIAL_CACHE_MB must be one of: ${LIBECPINT_RADIAL_CACHE_MB_VALUES}")
endif()

# configure the config header to pass the above variables to the program
configure_file (
  "${PROJECT_SOURCE_DIR}/include/libecpint/config.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/libecpint/config.hpp"
)

include(GNUInstallDirs)

if(MSVC)
    # MSVC does not include <cmath> constants, unless _USE_MATH_DEFINES is defined.
    add_definitions("/D_USE_MATH_DEFINES")
    # Set the exception handling model
    add_definitions("/EHsc")
endif()

add_subdirectory(src)

option(LIBECPINT_BUILD_TESTS "Enables Libecpint tests." ON)
if (LIBECPINT_BUILD_TESTS)
  message(STATUS "Will build libecpint tests.")
  include(CTest)
  include(external/ImportGTest.cmake)
  enable_testing()
  add_subdirectory(tests)
  add_subdirectory(benchmarks)
  add_subdirectory(example)
endif()

option(LIBECPINT_BUILD_DOCS "Enables Libecpint documentation." ON)
if (LIBECPINT_BUILD_DOCS)
  add_subdirectory(doc)
endif()

install(DIRECTORY include/libecpint DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} PATTERN "*.in" EXCLUDE)
install(FILES include/libecpint.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/include/libecpint/config.hpp" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libecpint)
install(DIRECTORY share/libecpint DESTINATION ${CMAKE_INSTALL_DATADIR})

message("\nSummary:")
message("--------")
message(STATUS "Maximum angular momentum unrolling : ${LIBECPINT_MAX_UNROL}")
message(STATUS "Radial factor cache per thread (MB): ${LIBECPINT_RADIAL_CACHE_MB}")
message(STATUS "Build tests                        : ${LIBECPINT_BUILD_TESTS}")
message("")
