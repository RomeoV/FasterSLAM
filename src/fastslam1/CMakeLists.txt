add_library(FasterSlamLib
  fastslam1.cpp
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
