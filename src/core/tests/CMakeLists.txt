# All test files get an executable and are registered as a test

file(GLOB testfiles "*_test.cpp")

foreach(file ${testfiles})
  get_filename_component(name ${file} NAME_WE)  # strip file ending
  # message("Detected test: ${name}")
  add_executable(core_${name}_exe
    ${file}
  )
  target_link_libraries(core_${name}_exe
    boost_ut
    FasterSlamCore
  )
  add_test(core_${name}
    core_${name}_exe
  )
endforeach()
