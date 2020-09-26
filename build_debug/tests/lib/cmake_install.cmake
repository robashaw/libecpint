# Install script for directory: /Users/robertshaw/devfiles/libecpint_new/tests/lib

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/robertshaw/devfiles/libecpint_new/build_debug/tests/lib/int_test1/cmake_install.cmake")
  include("/Users/robertshaw/devfiles/libecpint_new/build_debug/tests/lib/int_test2/cmake_install.cmake")
  include("/Users/robertshaw/devfiles/libecpint_new/build_debug/tests/lib/deriv_test1/cmake_install.cmake")
  include("/Users/robertshaw/devfiles/libecpint_new/build_debug/tests/lib/deriv_test2/cmake_install.cmake")
  include("/Users/robertshaw/devfiles/libecpint_new/build_debug/tests/lib/hess_test1/cmake_install.cmake")
  include("/Users/robertshaw/devfiles/libecpint_new/build_debug/tests/lib/hess_test2/cmake_install.cmake")
  include("/Users/robertshaw/devfiles/libecpint_new/build_debug/tests/lib/api_test1/cmake_install.cmake")
  include("/Users/robertshaw/devfiles/libecpint_new/build_debug/tests/lib/api_test2/cmake_install.cmake")
  include("/Users/robertshaw/devfiles/libecpint_new/build_debug/tests/lib/stress_test/cmake_install.cmake")

endif()

