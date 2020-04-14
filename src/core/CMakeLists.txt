add_library(FasterSlamCore
  example.cpp
  KF_cholesky_update.cpp
  multivariate_gauss.cpp
  linalg.cpp
  pi_to_pi.cpp
  predict_true.cpp
  stratified_random.cpp
  stratified_resample.cpp
  compute_steering.cpp
  add_observation_noise.cpp
  get_observations.cpp
)
target_include_directories(FasterSlamCore PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_subdirectory(tests)

if(TEST_COVERAGE)
  target_compile_options(FasterSlamCore PRIVATE
    --coverage)
  target_link_libraries(FasterSlamCore PRIVATE
    gcov)
endif()
