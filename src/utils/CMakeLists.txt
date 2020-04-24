add_library(rdtsc
  INTERFACE
)

target_include_directories(rdtsc
  INTERFACE 
  ./benchmarklib
)

add_library(FasterSlamUtils
  string_format.cpp
)

target_include_directories(FasterSlamUtils PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)


target_link_libraries(FasterSlamUtils
    rdtsc
  )

target_include_directories(FasterSlamUtils PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# add_subdirectory(tests)