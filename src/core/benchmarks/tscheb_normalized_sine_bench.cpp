#include "rdtsc_benchmark.h"
#include <iostream>
#include <cmath>
#include <functional>
#include "tscheb_sine.h"
#include "trigonometry.h"

#define NR 256

double sum = 0;
double* results;

template <typename FT>  // float type
void rand_angles_normalized(FT * m)
{
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<FT> dist(-M_PI, M_PI); // this is in [-PI, PI) so we negate later to get (-PI, PI]
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
    rand_angles_normalized(alphas);
}

void calc_normalized_csines(double* alphas) {
  double result = 0;
  for (size_t i = 0; i < NR; i++) {
      sum += sin(alphas[i]);
  }
}

void calc_tscheb_normalized_fsines(double* alphas) {
  float result = 0;
  for (size_t i = 0; i < NR; i++) {
      sum += tscheb_fsine(alphas[i], true);
  }
}

void calc_tscheb_normalized_dsines(double* alphas) {
  double result = 0;
  for (size_t i = 0; i < NR; i++) {
      sum += tscheb_dsine(alphas[i], true);
  }
}

void calc_vectorized_normalized_dsines(double* alphas) {
  tscheb_dsines(alphas, NR, results);
  sum += results[NR-1];
}

void calc_unrolled_normalized_dsines(double* alphas) {
  tscheb_dsines_unrolled(alphas, NR, results);
  sum += results[NR-1];
}

void calc_avx_normalized_dsines(double* alphas) {
  tscheb_dsines_avx(alphas, NR, results);
  sum += results[NR-1];
}

void read_sine(double* alphas) {
  float result = 0;
  for (size_t i = 0; i < NR; i++) {
      sum += read_sin(alphas[i]);
  }
}

void read_sine_vec(double* alphas) {
  float result = 0;
  for (size_t i = 0; i < NR; i+=4) {
      auto vec = _mm256_load_pd(alphas+i);
      _mm256_store_pd(results+i, read_sin2_vec(vec));
  }
  sum += results[NR-1];
}

int main() {
    // Initialize Input
    double *alphas_d;
    size_t N = NR;
    init_sin();
    init_sin2();

    build<double>(&alphas_d, N);

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&calc_tscheb_normalized_dsines)> bench("normalized_sine");

    // Set function to reload data for each benchmarked function
    bench.data_loader = &fill<double>;

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&calc_normalized_csines, "Cmath sines", 14*NR);
    bench.add_function(&calc_tscheb_normalized_fsines, "Tscheb. sines on norm. floats", 18*NR);
    bench.add_function(&calc_tscheb_normalized_dsines, "Tscheb. sines on norm. doubles", 18*NR);
    bench.add_function(&calc_vectorized_normalized_dsines, "Vect. tscheb. sines on norm. doubles", 18*NR);
    bench.add_function(&calc_unrolled_normalized_dsines, "Unrld. tscheb. sines on norm. doubles", 18*NR);
    bench.add_function(&calc_avx_normalized_dsines, "AVX. tscheb. sines on norm. doubles", 18*NR);
    bench.add_function(&read_sine, "read_sine", 6*NR);
    bench.add_function(&read_sine_vec, "read_sine_vec", 6*NR);

    // Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(alphas_d);

    // Output is given when bench is destroyed (to change this behaviour, set bench.destructor_output=false). Optionally you can call bench.summary() to get it.

    // If you want the summary to be written to a file, set bench.fout with your preferred ostream. (Default std::cout)

    // Free memory
    destroy(alphas_d);

    return 0;
}