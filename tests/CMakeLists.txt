# Disable a test (reported by ctest as "Not Run (Disabled)", not as a failure) when the library was
# built with LIBECPINT_MAX_L below what the test needs. Call after the test's add_test().
function(ecpint_test_requires_max_l test min_l)
  if (LIBECPINT_MAX_L LESS min_l)
    set_tests_properties(${test} PROPERTIES DISABLED TRUE)
    message(STATUS "Test ${test} disabled: needs LIBECPINT_MAX_L >= ${min_l} (have ${LIBECPINT_MAX_L})")
  endif()
endfunction()

add_subdirectory(units)
add_subdirectory(lib)

