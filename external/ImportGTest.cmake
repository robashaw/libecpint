find_library(GTEST_LIBRARY gtest
  HINTS "/usr/local" "/usr/local/lib" "/usr/lib"
)
find_path(GTEST_INCLUDE_DIR gtest/gtest.h
  HINTS "/usr/local/include" "/usr/include"
)

if((NOT GTEST_LIBRARY) OR (NOT GTEST_INCLUDE_DIR))
  message(STATUS "Unable to find google test, cloning and building ...")

  set(THREADS_REFER_PTHREAD_FLAG "ON")
  find_package(Threads REQUIRED)

  # Download and unpack googletest at configure time
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/external/CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()

  # Prevent overriding the parent project's compiler/linker
  # settings on Windows
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  # Set GTEST interface to temporarily installed libgtest.a
  set(GTEST_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/googletest-build/include)
  set(GTEST_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/googletest-build/lib/libgtest.a)
  add_library(gtest INTERFACE)
  target_include_directories(gtest INTERFACE ${GTEST_INCLUDE_DIR})
  target_link_libraries(gtest INTERFACE ${GTEST_LIBRARY} Threads::Threads)
else()
  message(STATUS "Found googletest: ${GTEST_INCLUDE_DIR} ${GTEST_LIBRARY}")
  add_library(gtest INTERFACE)
  target_include_directories(gtest INTERFACE ${GTEST_INCLUDE_DIR})
  target_link_libraries(gtest INTERFACE ${GTEST_LIBRARY})
endif()
