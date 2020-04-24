#include "estimate_json.h"
#include <iostream>
#include <array>
#include <algorithm>
#include <linalg.h>

#include "converters.h"
//C++ Code!

nlohmann::json estimate_step_json(Particle* particles, double* weights, int N, double time, int* ftag_visible, const size_t N_visible, int* da_table, int N_features, int id) {
    std::vector<double> weights_vec;
    std::vector<std::array<double,3>> particle_poses;

    std::vector<int> idf;
    std::vector<int> ftag_visible_vec;

    for (int i = 0; i< N_visible; i++) {
        ftag_visible_vec.push_back(ftag_visible[i]);
        idf.push_back(da_table[ftag_visible[i]]);
        //std::cout<<"idf["<<i<<"] = "<<da_table[ftag_visible[i]]<<" "<<ftag_visible[i]<<std::endl;
    }

    double x[3] = {0,0,0};
    double cov[9] = {0,0,0,0,0,0,0,0,0};
    

    double sum_w = 0.0;
    const int Nfa = particles[0].Nfa;
    double lms[Nfa*2];

    double lm_particles[N * 2*Nfa]; //row x col
    double lm_covs[2*Nfa * 2];

    for (int i = 0; i<Nfa; i++){
        lms[2*i] = 0;
        lms[2*i+1] = 0;

        lm_covs[4*i+0] = 0;
        lm_covs[4*i+1] = 0;
        lm_covs[4*i+2] = 0;
        lm_covs[4*i+3] = 0;

        for(int j = 0; j<N; j++) {
            lm_particles[2*Nfa*j + 2*i] = 0;
            lm_particles[2*Nfa*j + 2*i +1] = 0; 
        }
    }

    

    //print(lm_covs,2*Nfa,2,std::cout);

    for (int i = 0; i<N; i++) {
        double w = weights[i];
        weights_vec.push_back(w);
        x[0] += particles[i].xv[0] * w;
        x[1] += particles[i].xv[1] * w;
        x[2] += particles[i].xv[2] * w;
        particle_poses.push_back(Vector3dToStdVec(particles[i].xv));
        for (int j = 0; j < Nfa; j++) {
           
            lms[2*j] += particles[i].xf[2*j] * w;
            lms[2*j+1] += particles[i].xf[2*j+1] *w;

            lm_particles[2*Nfa*i + 2*j]       =   particles[i].xf[2*j]; //wrong here in indices but faster
            lm_particles[2*Nfa*i + 2*j + 1]   =   particles[i].xf[2*j+1]; 
        }
    }

    //print(lms,N_visible,2,std::cout);
    //print(lm_particles,N,2*Nfa,std::cout);


    
    for(int i = 0; i<N; i++){
        double w = weights[i];
        for (int j = 0; j < Nfa; j++) {
            
            double dx = lm_particles[2*Nfa*i + 2*j] - lms[2*j];
            double dy = lm_particles[2*Nfa*i + 2*j +1] - lms[2*j+1];

            lm_covs[4*j+0] += dx*dx*w;
            lm_covs[4*j+1] += dx*dy*w;
            lm_covs[4*j+2] += dy*dx*w;
            lm_covs[4*j+3] += dy*dy*w;
        }
    }
    //print(lms,N_visible,2,std::cout);
    std::vector<std::array<double,2>> visible_landmark_poses;
    std::vector<std::array<double,4>> visible_landmark_covs;

    for (int i = 0; i<N_visible;i++) {
        int index = idf[i];
        visible_landmark_poses.push_back(Vector2dToStdVec(lms+2*index));
        visible_landmark_covs.push_back(Matrix2dToStdVec(lm_covs+4*index));
    }

    return {
        {"id", id},
        {"timestamp",time},
        {"robot_pose",Vector3dToStdVec(x)},
        {"robot_covariance", Matrix3dToStdVec(cov)},
        {"particle_poses", particle_poses},
        {"weights", weights_vec},
        {"landmark_visible",ftag_visible_vec},
        {"landmark_poses", visible_landmark_poses},
        {"landmark_covariances", visible_landmark_covs}
    };
}   