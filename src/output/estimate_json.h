#pragma once
#include <iostream>
#include <array>
#include <vector>
#include "json.hpp"
#include "particle.h"


/*

namespace {
	auto Vector2dToStdVec = [](const Eigen::Vector2d& evec) -> std::array<double, 2> {
		return {evec(0), evec(1)};
	};
	auto Vector3dToStdVec = [](const Eigen::Vector3d& evec) -> std::array<double, 3> {
		return {evec(0), evec(1), evec(2)};
	};
	auto Matrix2dToStdVec = [](const Eigen::Matrix2d& evec) -> std::array<double, 4> {
		return {evec(0,0), evec(0,1), evec(1,0), evec(1,1)};
	};
    auto Matrix3dToStdVec = [](const Eigen::Matrix3d& evec) -> std::array<double, 9> {
		return {evec(0,0), evec(0,1), evec(0,2), evec(1,0), evec(1,1), evec(1,2),evec(2,0), evec(2,1), evec(2,2)};
	};
};

nlohmann::json estimate_step_json(const std::vector<Particle>& p_vec, double dt, std::vector<int>& ftag_visible, Eigen::VectorXd& da_table, int id);


*/