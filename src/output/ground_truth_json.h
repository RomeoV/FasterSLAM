#pragma once
#include <iostream>
#include <array>
#include <vector>
#include "json.hpp"


nlohmann::json ground_truth_step_json(Vector3d xtrue, double time, int id, int iwp, double G);

nlohmann::json ground_truth_keypoints_json(double* waypoints, double* landmarks, size_t N_w, size_t N_f);