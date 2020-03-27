project(FasterSLAM CXX)
cmake_minimum_required(VERSION 3.10)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE CACHE BOOL "Export compile commands" FORCE)

enable_testing()
add_subdirectory(external)
add_subdirectory(core)
# add_subdirectory(cpp/fastslam1)
