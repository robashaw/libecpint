
# ==================
#  Target: generate
# ==================

add_executable(generate
               generate.cpp
               lib/mathutil.cpp
               lib/angular.cpp)

target_compile_features(generate PUBLIC cxx_std_11)
target_include_directories(generate
  PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/libecpint
    ${CMAKE_CURRENT_BINARY_DIR}/../include/libecpint
)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/generated/ecpint_gen.part" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/generated/qgen.part" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/generated")
find_package (Python COMPONENTS Interpreter)
execute_process(COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/makelist.py" "${LIBECPINT_MAX_L}" "${CMAKE_CURRENT_BINARY_DIR}")
file(READ "${CMAKE_CURRENT_BINARY_DIR}/qlist.txt" GENERATED_SOURCES)
STRING(REGEX REPLACE "\n" ";" GENERATED_SOURCES "${GENERATED_SOURCES}")
add_custom_command(
	OUTPUT ${GENERATED_SOURCES}
	COMMAND generate "${PROJECT_SOURCE_DIR}/include/libecpint/"
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/generated/ecpint_gen.part
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

target_compile_features(ecpint PUBLIC cxx_std_11)
set_target_properties(ecpint PROPERTIES SOVERSION ${API_VERSION})
target_include_directories(ecpint
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/libecpint>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../include/libecpint>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/libecpint>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

if (LIBECPINT_USE_CERF)
    find_package(cerf CONFIG REQUIRED) # QUIET)  # COMPONENTS C)

    if(${cerf_FOUND})
        get_property(_loc TARGET cerf::cerf PROPERTY LOCATION)
        message(STATUS "${Cyan}Found cerf${ColourReset}: ${_loc}")
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_cerf})
            message(FATAL_ERROR "Suitable cerf could not be externally located as user insists")
        endif()
    endif()
    target_compile_definitions(ecpint PRIVATE USING_CERF=1)
    target_link_libraries(ecpint
     PRIVATE
      cerf::cerf
    )
else()
    message(STATUS "${Cyan}Using Faddeeva${ColourReset}: internal")
    target_link_libraries(ecpint
     PRIVATE
      Faddeeva
    )
endif()

if(LIBECPINT_USE_PUGIXML)
  target_include_directories(ecpint PRIVATE ${PUGIXML_INCLUDE_DIR})
  target_link_libraries(ecpint PRIVATE libpugixml)
endif()

# =========
#  Install
# =========

install(TARGETS ecpint EXPORT ecpintTargets DESTINATION ${CMAKE_INSTALL_LIBDIR})
if (LIBECPINT_USE_PUGIXML)
  # used in api.h
  target_compile_definitions(ecpint PUBLIC HAS_PUGIXML)
endif()
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

install(TARGETS ecpint DESTINATION ${CMAKE_INSTALL_LIBDIR})
