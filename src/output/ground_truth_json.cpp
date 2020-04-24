#include "estimate_json.h" 
#include "ground_truth_json.h" 
#include <iostream>

/*

nlohmann::json ground_truth_step_json(Eigen::Vector3d xtrue, double time, int id, int iwp, double G) {

    return {
        {"id", id},
        {"timestamp",time},
        {"robot_pose",Vector3dToStdVec(xtrue)},
        {"current_waypoint",iwp},
        {"Steering angle", G}
    };
}   



nlohmann::json ground_truth_keypoints_json(Eigen::MatrixXd& waypoints, Eigen::MatrixXd& landmarks) {
    std::vector<std::array<double,2>> landmark_poses;
    std::vector<std::array<double,2>> waypoint_poses;

    for (int i = 0; i<waypoints.cols();i++) {
        waypoint_poses.push_back(Vector2dToStdVec(waypoints.col(i)));
    }


    for (int i = 0; i<landmarks.cols();i++) {
        landmark_poses.push_back(Vector2dToStdVec(landmarks.col(i)));
    }

    return {
        {"waypoints", waypoint_poses},
        {"landmarks", landmark_poses}
    };

}

*/