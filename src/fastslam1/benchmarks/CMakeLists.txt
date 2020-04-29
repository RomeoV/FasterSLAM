# All bench files get an executable and are registered as a bench

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
  FasterSlamLib
)

# add_test(core_${name}
# fastslam1_${name}_exe
# )

endforeach()