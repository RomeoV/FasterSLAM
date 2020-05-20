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

    // Initialize Input
    double *alphas_d;
    size_t N = NR;

    build<double>(&alphas_d, N);

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&calc_tscheb_dsines)> bench("not_normalized_sine");

    // Set function to reload data for each benchmarked function
    bench.data_loader = &fill<double>;

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&calc_csines, "Cmath sines", 14*NR);
    bench.add_function(&calc_tscheb_fsines, "Tscheb. sines on floats", (18+10)*NR);  // the 10 is only approximate and should probably be performance counted
    bench.add_function(&calc_tscheb_dsines, "Tscheb. sines on doubles", (18+10)*NR);
    bench.add_function(&calc_vectorized_dsines, "Vect. tscheb. sines on doubles", 18*NR);
    bench.add_function(&calc_unrolled_dsines, "Unrld. tscheb. sines on doubles", 18*NR);
    bench.add_function(&calc_avx_dsines, "AVX tscheb. sines on doubles", 18*NR);

    // Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(alphas_d);

    // Output is given when bench is destroyed (to change this behaviour, set bench.destructor_output=false). Optionally you can call bench.summary() to get it.

    // If you want the summary to be written to a file, set bench.fout with your preferred ostream. (Default std::cout)

    // Free memory
    destroy(alphas_d);

    return 0;
}
