add_library(FasterSlamCore
  KF_cholesky_update.cpp
  add_feature.cpp
  add_observation_noise.cpp
  compute_jacobians.cpp
  compute_steering.cpp
  compute_weight.cpp
  configfile.cpp
  example.cpp
  tscheb_sine.cpp
  feature_update.cpp
  get_observations.cpp
  linalg.cpp
  multivariate_gauss.cpp
  particle.cpp
  pi_to_pi.cpp
  predict.cpp
  predict_true.cpp
  resample_particles.cpp
  stratified_random.cpp
  stratified_resample.cpp
  add_control_noise.cpp
  data_associate_known.cpp
  fastrand.cpp
  trigonometry.cpp
)
target_include_directories(FasterSlamCore PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_options(FasterSlamCore PUBLIC
  -march=native
)

add_subdirectory(tests)
add_subdirectory(benchmarks)

if(TEST_COVERAGE)
  target_compile_options(FasterSlamCore PRIVATE
    --coverage)
  target_link_libraries(FasterSlamCore PRIVATE
    gcov)
endif()
