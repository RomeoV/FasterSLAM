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
  compute_weight_and_feature_update.cpp
)
if (USE_YGLEE_KF_CHOLESKY)
    target_compile_options(FasterSlamCore PUBLIC -DKF_YGLEE)
endif()
target_include_directories(FasterSlamCore PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

if(${USE_MARCH_NATIVE})
target_compile_options(FasterSlamCore
  PUBLIC
  -march=native
  -O3
  -ffast-math
)
endif()

target_compile_options(FasterSlamCore
  PUBLIC
  -Wall
  -Wno-sign-compare
)

add_subdirectory(tests)
if(${BUILD_BENCHMARKS})
add_subdirectory(benchmarks)
endif()

if(TEST_COVERAGE)
  target_compile_options(FasterSlamCore PRIVATE
    --coverage)
  target_link_libraries(FasterSlamCore PRIVATE
    gcov)
endif()
