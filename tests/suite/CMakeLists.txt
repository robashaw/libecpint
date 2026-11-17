# Integral checks over the systems of tests/data/systems.txt, one ctest entry per (check, system):
# Suite.Stored.<system>, Suite.Screening.<system>, Suite.Threads.<system>, and with
# LIBECPINT_TEST_REFERENCE=ON Suite.Live.<system> (label "reference") for systems with live=yes.
option(LIBECPINT_TEST_REFERENCE
       "Add the (slow) tests against the reference engine computed at test time, label 'reference'"
       OFF)

add_executable(suite_test suite.cpp)
target_link_libraries(suite_test ecpint_testsupport)

# Regenerates every stored reference file from the reference engine
add_executable(regenerate_refs regenerate.cpp)
target_link_libraries(regenerate_refs ecpint_testsupport)
add_custom_target(regenerate_references
  COMMAND regenerate_refs
  DEPENDS regenerate_refs
  USES_TERMINAL
  COMMENT "Writing tests/data/reference from the reference engine")

set(MANIFEST ${PROJECT_SOURCE_DIR}/tests/data/systems.txt)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${MANIFEST})
file(STRINGS ${MANIFEST} manifest_lines)
foreach(line IN LISTS manifest_lines)
  string(REGEX REPLACE "#.*" "" line "${line}")
  if(NOT line MATCHES "name=([^ \t]+)")
    continue()
  endif()
  set(name ${CMAKE_MATCH_1})
  set(checks Stored Screening Threads)
  if(LIBECPINT_TEST_REFERENCE AND line MATCHES "live=yes")
    list(APPEND checks Live)
  endif()
  foreach(check IN LISTS checks)
    add_test(NAME Suite.${check}.${name}
             COMMAND suite_test --gtest_filter=Systems/${check}.*/${name})
    if(check STREQUAL "Live")
      set_tests_properties(Suite.${check}.${name} PROPERTIES LABELS reference)
    endif()
  endforeach()
endforeach()
