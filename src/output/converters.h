#pragma once

#include <vector>
#include <array>


namespace {
	auto Vector2dToStdVec = [](double* evec) -> std::array<double, 2> {
		return {evec[0], evec[1]};
	};
	auto Vector3dToStdVec = [](double* evec) -> std::array<double, 3> {
		return {evec[0], evec[1], evec[2]};
	};
	auto Matrix2dToStdVec = [](double* evec) -> std::array<double, 4> {
		return {evec[0], evec[1], evec[2], evec[3]};
	};
    auto Matrix3dToStdVec = [](double* evec) -> std::array<double, 9> {
		return {evec[0], evec[1], evec[2], evec[3], evec[4], evec[5],evec[6], evec[7], evec[8]};
	};
};