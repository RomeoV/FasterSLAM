#include "estimate_json.h" 
#include "ground_truth_json.h" 
#include <iostream>

#include "converters.h"

nlohmann::json ground_truth_step_json(Vector3d xtrue, double time, int id, int iwp, double G) {
    return {
        {"id", id},
        {"timestamp",time},
        {"robot_pose",Vector3dToStdVec(xtrue)},
        {"current_waypoint",iwp},
        {"Steering angle", G}
    };
}   



nlohmann::json ground_truth_keypoints_json(double* waypoints, double* landmarks, size_t N_w, size_t N_f) {
    std::vector<std::array<double,2>> landmark_poses;
    std::vector<std::array<double,2>> waypoint_poses;



    for (int i = 0; i<N_w;i++) {
        waypoint_poses.push_back(Vector2dToStdVec(waypoints+2*i));
    }


    for (int i = 0; i<N_f;i++) {
        landmark_poses.push_back(Vector2dToStdVec(landmarks+2*i));
    }

    return {
        {"waypoints", waypoint_poses},
        {"landmarks", landmark_poses}
    };

}