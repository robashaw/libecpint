# Shared code for the integral tests: loading test systems, the independent reference engine with
# the library's assembly, comparisons, and stored reference files.
find_package(Threads REQUIRED)

execute_process(
  COMMAND git describe --always --dirty
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE LIBECPINT_GIT_DESCRIBE
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NOT LIBECPINT_GIT_DESCRIBE)
  set(LIBECPINT_GIT_DESCRIBE "unknown")
endif()

add_library(ecpint_testsupport STATIC
  system.cpp
  reference.cpp
  compare.cpp
  reference_io.cpp)
target_include_directories(ecpint_testsupport PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ecpint_testsupport PUBLIC ecpint gtest Threads::Threads)
target_compile_features(ecpint_testsupport PUBLIC cxx_std_17)
target_compile_definitions(ecpint_testsupport PUBLIC
  LIBECPINT_TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/tests/data"
  LIBECPINT_SHARE_DIR="${PROJECT_SOURCE_DIR}/share/libecpint"
  LIBECPINT_TEST_GENERATOR="libecpint-${LIBECPINT_GIT_DESCRIBE}")
