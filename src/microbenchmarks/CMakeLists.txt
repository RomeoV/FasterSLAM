file(GLOB benchfiles "*_bench.cpp")

if(NOT ${CMAKE_BUILD_TYPE} STREQUAL Release)
    set(bench_comment "Please make sure to set build type to Release for benchmarking!")
endif()
add_custom_command(OUTPUT _benchmark_command
    COMMAND echo "Starting benchmarks"
    COMMENT ${bench_comment}
)

foreach(file ${benchfiles})
  get_filename_component(name ${file} NAME_WE)  # strip file ending
  # message("Detected benchmark: ${name}")
  add_executable(${name}_exe
    ${file}
  )
  target_link_libraries(${name}_exe
    FasterSlamCore
    boost_ut
    nanobench
  )
  add_custom_command(OUTPUT _benchmark_command
      COMMAND ${name}_exe
      APPEND
  )
endforeach()

add_custom_target(benchmarks
    DEPENDS _benchmark_command
)
