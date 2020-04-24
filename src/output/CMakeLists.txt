add_library(FasterSlamOutput
  estimate_json.cpp
  ground_truth_json.cpp
  converters.h
)

target_include_directories(FasterSlamOutput PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)


target_link_libraries(FasterSlamOutput
    nlohmann_json
    FasterSlamCore
  )

target_include_directories(FasterSlamOutput PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_subdirectory(tests)

add_executable(FasterSlamVisualize
  visualize.cpp
)
target_link_libraries(FasterSlamVisualize PUBLIC
  FasterSlamCore
  FasterSlamLib
  nlohmann_json
  FasterSlamOutput
)