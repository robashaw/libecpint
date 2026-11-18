
# ==================
#  Target: generate
# ==================

if(CMAKE_CROSSCOMPILING)
    find_package(genecpint)
else()
    add_executable(generate
                   generate.cpp
                   lib/mathutil.cpp
                   lib/angular.cpp)

    target_compile_features(generate PUBLIC cxx_std_17)
    target_include_directories(generate
      PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/include/libecpint
        ${CMAKE_CURRENT_BINARY_DIR}/../include/libecpint
      )
    export(
      TARGETS
        generate
      FILE
        "${CMAKE_BINARY_DIR}/genecpintConfig.cmake"
      )
endif()

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
find_package (Python COMPONENTS Interpreter)
execute_process(COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/makelist.py" "${LIBECPINT_MAX_UNROL}" "${CMAKE_CURRENT_BINARY_DIR}")
file(READ "${CMAKE_CURRENT_BINARY_DIR}/qlist.txt" GENERATED_SOURCES)
STRING(REGEX REPLACE "\n" ";" GENERATED_SOURCES "${GENERATED_SOURCES}")
# The generated header qgen.hpp is internal, and written to the build tree with the sources
add_custom_command(
	OUTPUT ${GENERATED_SOURCES} "${CMAKE_CURRENT_BINARY_DIR}/generated/qgen.hpp"
	COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/generated/ecpint_gen.part"
	        "${CMAKE_CURRENT_SOURCE_DIR}/generated/qgen.part" "${CMAKE_CURRENT_BINARY_DIR}/generated"
	COMMAND generate "${CMAKE_CURRENT_BINARY_DIR}/generated/"
	DEPENDS generate ${CMAKE_CURRENT_SOURCE_DIR}/generated/ecpint_gen.part
	        ${CMAKE_CURRENT_SOURCE_DIR}/generated/qgen.part
	COMMENT "Generating and unrolling integral code."
)

# ================
#  Target: ecpint
# ================

file(GLOB_RECURSE LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cpp")
add_library(ecpint ${LIB_SOURCES} ${GENERATED_SOURCES})
add_library(ECPINT::ecpint ALIAS ecpint)

if(NOT BUILD_SHARED_LIBS)
  set_property(TARGET ecpint PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

target_compile_features(ecpint PUBLIC cxx_std_17)
set_target_properties(ecpint PROPERTIES SOVERSION ${API_VERSION})
target_include_directories(ecpint
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/libecpint>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../include/libecpint>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/libecpint>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
# BEFORE, so a stale qgen.hpp left in include/libecpint by the old generator is never picked up
target_include_directories(ecpint BEFORE PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)

# =========
#  Install
# =========

install(
  TARGETS
    ecpint
  EXPORT
    ecpintTargets
  RUNTIME DESTINATION
    ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION
    ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION
    ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}
  )
include(CMakePackageConfigHelpers)
# Version file
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/ecpint-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)
# Config file
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/ecpint-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ecpint
)
# Targets files
export(
  EXPORT ecpintTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/ecpint-targets.cmake
)
install(
  EXPORT ecpintTargets
  FILE ecpint-targets.cmake
  NAMESPACE ECPINT::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ecpint
)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/ecpint-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/ecpint-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ecpint
)

