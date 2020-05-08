# All bench files get an executable and are registered as a bench

add_custom_command(OUTPUT _run_benchmarks
    POST_BUILD
    COMMAND echo "Starting benchmarks"
)

file(GLOB benchfiles "*_bench.cpp")
  foreach(file ${benchfiles})
    get_filename_component(name ${file} NAME_WE)  # strip file ending
    # message("Detected Benchmark: ${name}")
    add_executable(core_${name}_exe
      ${file}
    )
  target_link_libraries(core_${name}_exe
    boost_ut
    FasterSlamUtils
    FasterSlamCore
  )

  add_custom_command(OUTPUT _run_benchmarks
    COMMAND core_${name}_exe
    APPEND
  )
endforeach()

add_custom_target(benchmarks DEPENDS _run_benchmarks)