#pragma once

#include <math.h>
#include <iostream>
#include <immintrin.h>
const int n_2pi = 5;
const int steps_pi = 1000;

const double step_size = M_PI / steps_pi;
const __m256d step_size_vec = _mm256_set1_pd(step_size);
const __m256d minus_ones = _mm256_set1_pd(-1.0);
const __m256d pi_2_vec = _mm256_set1_pd(M_PI_2);

const __m256d offset_sin2 = _mm256_set1_pd(n_2pi * M_PI * step_size);

const int N_sin = n_2pi * (steps_pi*2 +1);

// One-sided sine table [- n_2pi * M_PI, + n_2pi * M_PI]
double sin_table[N_sin]; //Here we only store positive sin for angles >=0


// Symmetric sine table [- n_2pi * M_PI, + n_2pi * M_PI]
double sin2_table[N_sin];

void init_sin() {
    for (int i = 0; i < N_sin; i++) {
        sin_table[i] = std::sin((i/steps_pi)*M_PI);
    }
}

void init_sin2() {
    //symmetric
    double low = - n_2pi * M_PI;
    for (int i = 0; i < N_sin; i++) {
        sin2_table[i] = std::sin(low + (i/steps_pi)*M_PI);
    }
}



double read_sin(double angle) {
    int idx = int(angle / M_PI * steps_pi);
    if (angle >= 0) {
        return sin_table[idx];
    } else {
        return -sin_table[-idx];
    }
}

inline __m256d read_sin_vec(__m256d angle) {
    angle = _mm256_mul_pd(step_size_vec, angle);
    __m128i index = _mm256_cvtpd_epi32(angle);

    //Alternative:
    //__m128i all_pos_index = _mm_sign_epi32(index, index);
    //auto load_indices64 = _mm256_cvtepi32_epi64(all_pos_index);

    auto load_indices64 = _mm256_cvtepu32_epi64(index);

    auto results = _mm256_maskload_pd(sin_table, load_indices64);
    auto neg_results = _mm256_mul_pd(minus_ones, results);
    return _mm256_blendv_pd(results, neg_results,angle);
}

inline __m256d read_cos_vec(__m256d angle) {
    //This sub and the mul in read_sin_vec can be fused to an fmsub, but I keep it modular for now
    // In cos we can also use symmetry!
    angle = _mm256_sub_pd(angle, pi_2_vec);
    return read_sin_vec(angle);
}

inline __m256d read_sin2_vec(__m256d angle) {
    angle = _mm256_fmadd_pd(step_size_vec, angle, offset_sin2);
    auto load_indices64 = _mm256_castpd_si256(angle);
    return _mm256_maskload_pd(sin2_table, load_indices64);
}

inline __m256d read_cos2_vec(__m256d angle) {
    angle = _mm256_sub_pd(angle, pi_2_vec);
    return read_sin2_vec(angle);
}

double read_cos(double angle) {
    return read_sin(angle - M_PI_2);
}

