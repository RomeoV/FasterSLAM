add_library(boost_ut
  INTERFACE
)
target_include_directories(boost_ut
  INTERFACE
  ./ut
)
target_compile_features(boost_ut
  INTERFACE
  cxx_std_17
)

add_library(nanobench
  nanobench/nanobench.cpp
)
target_include_directories(nanobench
  PUBLIC
  ./nanobench
)
target_compile_features(nanobench
  PUBLIC
  cxx_std_17
)

add_library(nlohmann_json
  INTERFACE
)

target_include_directories(nlohmann_json 
  INTERFACE 
  ./json
)
target_compile_features(nlohmann_json
  INTERFACE
  cxx_std_17
)