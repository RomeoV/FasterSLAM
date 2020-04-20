#include "estimate_json.h"
#include <iostream>

/*

nlohmann::json estimate_step_json(const std::vector<Particle>& p_vec, double time, std::vector<int>& ftag_visible, Eigen::VectorXd& da_table, int id) {
    std::vector<double> weights;
    std::vector<std::array<double,3>> particle_poses;

    std::vector<int> idf;

    for (auto ftag : ftag_visible) {
        idf.push_back(da_table(ftag));
    }
    

    double sum_w = 0.0;
    Eigen::Vector3d x;
    Eigen::Matrix3d cov;
    const int num_landmarks = static_cast<int>(p_vec.at(0).xf().size());;
    Eigen::MatrixXd lms, lm_covs;
    lms.resize(2,num_landmarks);
    lm_covs.resize(2,2*num_landmarks);
    x.setZero();
    cov.setZero();
    lms.setZero();
    lm_covs.setZero();
    std::vector<double> w;

    for (auto part : p_vec) {
        sum_w+=part.w();
        x=x+part.xv() * part.w();
        cov = cov + part.Pv() * part.w();
        weights.push_back(part.w());
        particle_poses.push_back(Vector3dToStdVec(part.xv()));
        for (int j = 0; j<num_landmarks;j++) {
            lms.col(j) += part.xf().at(j) * part.w();
            lm_covs.block<2,2>(0,2*j)+=part.Pf().at(j) * part.w();
        }
    }
    lms = lms/sum_w;
    lm_covs = lm_covs / sum_w;
    std::vector<std::array<double,2>> visible_landmark_poses;
    std::vector<std::array<double,4>> visible_landmark_covs;

    for (int i = 0; i<ftag_visible.size();i++) {
        visible_landmark_poses.push_back(Vector2dToStdVec(lms.col(idf.at(i))));
        visible_landmark_covs.push_back(Matrix2dToStdVec(lm_covs.block<2,2>(0,2*idf.at(i))));
    }

    return {
        {"id", id},
        {"timestamp",time},
        {"robot_pose",Vector3dToStdVec(x)},
        {"robot_covariance", Matrix3dToStdVec(cov)},
        {"particle_poses", particle_poses},
        {"weights", weights},
        {"landmark_visible",ftag_visible},
        {"landmark_poses", visible_landmark_poses},
        {"landmark_covariances", visible_landmark_covs}
    };
}   

*/