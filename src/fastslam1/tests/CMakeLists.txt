# All test files get an executable and are registered as a test

file(GLOB testfiles "*_test.cpp")

foreach(file ${testfiles})
  get_filename_component(name ${file} NAME_WE)  # strip file ending
  # message("Detected test: ${name}")
  add_executable(sim_${name}_exe
    ${file}
  )
  target_link_libraries(sim_${name}_exe
    boost_ut
    FasterSlamCore
    FasterSlamLib
  )
  add_test(sim_${name}
    sim_${name}_exe
  )
endforeach()

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/inputfiles_test)
configure_file(
    test_webmap.mat
    ${CMAKE_CURRENT_BINARY_DIR}/inputfiles_test/test_webmap.mat
    COPYONLY
)
