add_library(FasterSlamCore
  add_feature.cpp
  compute_jacobians.cpp
  compute_steering.cpp
  compute_weight.cpp
  example.cpp
  feature_update.cpp
  KF_cholesky_update.cpp
  linalg.cpp
  multivariate_gauss.cpp
  particle.cpp
  pi_to_pi.cpp
  predict.cpp
  predict_true.cpp
  resample_particles.cpp
  stratified_random.cpp
  stratified_resample.cpp
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
