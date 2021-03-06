# All test files get an executable and are registered as a test

file(GLOB testfiles "*_test.cpp")

foreach(file ${testfiles})
  get_filename_component(name ${file} NAME_WE)  # strip file ending
  # message("Detected test: ${name}")
  add_executable(core_${name}_exe
    ${file}
  )
  if (USE_YGLEE_KF_CHOLESKY)
      target_compile_options(core_${name}_exe PUBLIC -DKF_YGLEE)
  endif()
  target_link_libraries(core_${name}_exe
    boost_ut
    FasterSlamCore
    FasterSlamUtils
  )
  add_test(core_${name}
    core_${name}_exe
  )
endforeach()

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/inputfiles_test)
configure_file(
    test_webmap.mat
    ${CMAKE_CURRENT_BINARY_DIR}/inputfiles_test/test_webmap.mat
    COPYONLY
)
configure_file(
    lm.txt
    ${CMAKE_CURRENT_BINARY_DIR}/inputfiles_test/lm.txt
    COPYONLY
)
