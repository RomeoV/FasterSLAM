add_library(FasterSlamCore
  example.cpp
)
target_include_directories(FasterSlamCore PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

if(TEST_COVERAGE)
  include(CodeCoverage)
  setup_target_for_coverage_gcovr_html(
    NAME core_coverage
    EXECUTABLE ctest .
    DEPENDENCIES FasterSlamCore 
    EXCLUDE ${PROJECT_SOURCE_DIR}/external ${PROJECT_SOURCE_DIR}/core/tests
    BASE_DIRECTORY "${PROJECT_SOURCE_DIR}/core"
  )
  target_link_libraries(FasterSlamCore gcov)
  append_coverage_compiler_flags()
endif()

add_subdirectory(tests)
