#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "configfile.h"
#include "compute_steering.h"
#include "predict_true.h"
#include "particle.h"

/*!
    Calculates the particles and their positions. Mostly calls other functions.
    @param[out] particles  All particles
    @param[in]  lm        list of landmark data
    @param[in]  wp        list of waypoints, only used to compute steering.
 */
void fastslam1_sim(double* lm, const size_t lm_rows, const size_t lm_cols, 
                   double* wp, const size_t wp_rows, const size_t wp_cols, 
                   Particle **particles_, double **weights_);

void fastslam1_sim_active( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_);

void fastslam1_sim_base( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_);

double fastslam1_sim_base_flops( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_);

double fastslam1_sim_base_memory( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_);

double fastslam1_sim_active_flops( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_);

double fastslam1_sim_active_memory( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_);
    
