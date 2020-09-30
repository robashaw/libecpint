add_library(Faddeeva Faddeeva.cpp)
set_property(TARGET Faddeeva PROPERTY CXX_STANDARD 11)
target_include_directories(Faddeeva
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/libecpint>
)
install(TARGETS Faddeeva EXPORT ecpintTargets DESTINATION lib)
