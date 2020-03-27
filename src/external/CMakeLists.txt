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
# set_target_properties(boost_ut
#   PROPERTIES
#   CXX_STANDARD 17
# )
