#pragma once

#include <math.h>
#include <iostream>
#include <immintrin.h>


#define n_2pi 3
#define steps_pi 10000
#define N_sin n_2pi * (steps_pi*2 +1)

extern const double step_size;

extern const double num_steps;

extern const __m256d step_size_vec;
extern const __m256d num_steps_vec;
extern const __m256d minus_ones;
extern const __m256d pi_2_vec;

extern const __m256d offset_sin2;
// One-sided sine table [- n_2pi * M_PI, + n_2pi * M_PI]
extern double sin_table[N_sin]; //Here we only store positive sin for angles >=0


// Symmetric sine table [- n_2pi * M_PI, + n_2pi * M_PI]
extern double sin2_table[N_sin];





void init_sin();

void init_sin2();


double read_sin(double angle);

//! ERROR in 0.x range (rad)!
float atan2_approximation1(float y, float x);

//! Error in 0.00x range (appr)
float sqrt3(const float& n);

float rsqrtss_times_x( float* x);
void rsqrtss_times_x_void( float* x, float* out);

__m256d read_sin_vec(__m256d angle);
__m256d read_cos_vec(__m256d angle);
__m256d read_sin2_vec(__m256d angle);
__m256d read_cos2_vec(__m256d angle);
double read_cos(double angle);
double Q_rsqrt( double number );