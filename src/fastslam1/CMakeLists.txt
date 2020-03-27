add_library(FasterSlamLib

)
target_include_directories(FasterSlamLib PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(FasterSlamLib PUBLIC
  FasterSlamCore
)

add_subdirectory(tests)
