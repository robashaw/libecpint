add_library(Faddeeva Faddeeva.cpp)
set_property(TARGET Faddeeva PROPERTY CXX_STANDARD 11)
target_include_directories(Faddeeva
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/libecpint>
)
if(NOT BUILD_SHARED_LIBS)
  set_property(TARGET Faddeeva PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
install(TARGETS Faddeeva EXPORT ecpintTargets DESTINATION ${CMAKE_INSTALL_LIBDIR})
