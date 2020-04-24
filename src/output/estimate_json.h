#pragma once
#include <iostream>
#include <array>
#include <vector>
#include "json.hpp"
#include "particle.h"


nlohmann::json estimate_step_json(Particle* particles, double* weights, int N, double time, int* ftag_visible, size_t N_visible, int* da_table, int N_features, int id)  ;