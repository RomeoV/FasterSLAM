#include "rdtsc_benchmark.h"
#include <iostream>
#include <cmath>
#include <functional>
#include "tscheb_sine.h"
#include "pi_to_pi.h"

#define NR 256

double sum = 0;
double* results;

template <typename FT>  // float type
void rand_angles(FT * m)
{
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<FT> dist(-10*M_PI, 10*M_PI); // this is in [-PI, PI) so we negate later to get (-PI, PI]
    for (size_t i = 0; i < NR; ++i)  
        m[i] = -dist(gen);
}

template <typename FT>
void build(FT **a, size_t N)
{
    *a = static_cast<FT *>(aligned_alloc(32, N * sizeof(FT)));
    results = static_cast<double*>(aligned_alloc(32, N*sizeof(double)));
}

template <typename FT>
void destroy(FT * m)
{
    free(m);
    free(results);
}

template <typename FT>
void fill (FT* alphas) {
    rand_angles(alphas);
}

void calc_csines(double* alphas) {
  for (size_t i = 0; i < NR; i++) {
      sum += sin(alphas[i]);
  }
}

void calc_tscheb_fsines(double* alphas) {
  for (size_t i = 0; i < NR; i++) {
      sum += tscheb_fsine(alphas[i], false);
  }
}

void calc_tscheb_dsines(double* alphas) {
  for (size_t i = 0; i < NR; i++) {
      sum += tscheb_dsine(alphas[i], false);
  }
}

void calc_vectorized_dsines(double* alphas) {
  double normalized_alphas[NR];
  for (size_t i = 0; i < NR; i++) {
      normalized_alphas[i] = pi_to_pi_while(alphas[i]);
  }
  tscheb_dsines(normalized_alphas, NR, results);
  sum += results[NR-1];
}

void calc_unrolled_dsines(double* alphas) {
  double normalized_alphas[NR];
  for (size_t i = 0; i < NR; i++) {
      normalized_alphas[i] = pi_to_pi_while(alphas[i]);
  }
  tscheb_dsines_unrolled(normalized_alphas, NR, results);
  sum += results[NR-1];
}

void calc_avx_dsines(double* alphas) {
  double normalized_alphas[NR];
  for (size_t i = 0; i < NR; i++) {
      normalized_alphas[i] = pi_to_pi_while(alphas[i]);
  }
  tscheb_dsines_avx(normalized_alphas, NR, results);
  sum += results[NR-1];
}

int main() {

    // not used
    // seems mostly to be the same as tsched_normalized_sine_bench

    return 0;
}
