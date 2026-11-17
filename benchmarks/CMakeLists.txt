# Timing benchmark, not run by ctest (see README.md): build and run with the Benchmark target, or
# compare results between builds with compare.py
add_executable(ecpint_bench bench.cpp)
target_link_libraries(ecpint_bench ecpint_testsupport)

add_custom_target(Benchmark
  COMMAND ecpint_bench
  DEPENDS ecpint_bench
  USES_TERMINAL)
