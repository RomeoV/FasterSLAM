project(FasterSLAM CXX)
cmake_minimum_required(VERSION 3.10)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE CACHE BOOL "Export compile commands" FORCE)

set(TEST_COVERAGE FALSE CACHE BOOL "Debug compile and test coverage output")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

include(CTest)
enable_testing()
add_subdirectory(core)
add_subdirectory(fastslam1)
add_subdirectory(external)
add_subdirectory(microbenchmarks)
add_subdirectory(output)
add_subdirectory(utils)

if(TEST_COVERAGE)
  include(TestCoverage)
endif()
