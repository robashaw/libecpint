if (NOT LIBECPINT_USE_CERF)
add_library(Faddeeva Faddeeva.cpp)
target_compile_features(Faddeeva PUBLIC cxx_std_11)
set_target_properties(Faddeeva PROPERTIES SOVERSION ${API_VERSION})
target_include_directories(Faddeeva
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/libecpint>
)
if(NOT BUILD_SHARED_LIBS)
  set_property(TARGET Faddeeva PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
install(
  TARGETS
    Faddeeva
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
endif()
