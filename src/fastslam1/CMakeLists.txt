add_library(FasterSlamLib
  fastslam1_utils.cpp
  predict_update.cpp
  observe_update.cpp
  fastslam1_sim.cpp
)
target_include_directories(FasterSlamLib PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(FasterSlamLib PUBLIC
  FasterSlamCore
  FasterSlamUtils
)

add_executable(FasterSlamExe
  main.cpp
)
target_link_libraries(FasterSlamExe PUBLIC
  FasterSlamCore
  FasterSlamLib
)

add_subdirectory(tests)
add_subdirectory(benchmarks)

target_compile_options(FasterSlamCore PRIVATE -march=core-avx-i)
target_compile_options(FasterSlamCore PRIVATE -march=core-avx2)
target_compile_options(FasterSlamCore PRIVATE -march=skylake-avx512)

if(TEST_COVERAGE)
  target_compile_options(FasterSlamLib PRIVATE
    --coverage
  )
  target_link_libraries(FasterSlamLib PRIVATE
    gcov
  )
endif()
