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
)

add_executable(FasterSlamExe
  main.cpp
)
target_link_libraries(FasterSlamExe PUBLIC
  FasterSlamCore
  FasterSlamLib
)

add_subdirectory(tests)

if(TEST_COVERAGE)
  target_compile_options(FasterSlamLib PRIVATE
    --coverage
  )
  target_link_libraries(FasterSlamLib PRIVATE
    gcov
  )
endif()
