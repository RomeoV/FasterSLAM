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

if(TEST_COVERAGE)
  include(CodeCoverage)
  setup_target_for_coverage_gcovr_html(
    NAME lib_coverage
    BASE_DIRECTORY "${PROJECT_SOURCE_DIR}/fastslam1"
    EXECUTABLE ctest .
    DEPENDENCIES FasterSlamLib 
    EXCLUDE ${PROJECT_SOURCE_DIR}/external ${PROJECT_SOURCE_DIR}/fastslam1/tests
  )
  target_link_libraries(FasterSlamLib PUBLIC gcov)
  append_coverage_compiler_flags()
endif()

