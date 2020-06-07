#include "rdtsc_benchmark.h"
#include <iostream>
#include <cmath>
#include <functional>
#include "typedefs.h"
#include "tscheb_sine.h"
#include "trigonometry.h"
#include "pi_to_pi.h"

#include "ut.hpp"
using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

size_t NR = 256;
size_t NR_max = pow(2,16);

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
  for (size_t i = 0; i < NR; i++) {
      sum += sin(alphas[i]);
  }
}

void calc_tscheb_normalized_fsines(double* alphas) {
  for (size_t i = 0; i < NR; i++) {
      sum += tscheb_fsine(alphas[i], true);
  }
}

void calc_tscheb_normalized_dsines(double* alphas) {
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
  for (size_t i = 0; i < NR; i++) {
      sum += read_sin(alphas[i]);
  }
}

#ifdef __AVX2__
void read_sine_vec(double* alphas) {
  for (size_t i = 0; i < NR; i+=4) {
      auto vec = _mm256_load_pd(alphas+i);
      _mm256_store_pd(results+i, read_sin2_vec(vec));
  }
  sum += results[NR-1];
}
#endif

int main() {
    // Initialize Input
    double *alphas_d;
    size_t N = NR;
    init_sin();
    init_sin2();

    build<double>(&alphas_d, NR_max);

    // sanity checks
    double sin_results[99];
    double fsine_results[99];
    double dsine_results[99];
    double dsine_vec_results[99];
    double dsine_unrolled_results[99];
    double* angles = (double*)aligned_alloc(32, 100*sizeof(double));
    double* angles_normalized = (double*)aligned_alloc(32, 100*sizeof(double));
    auto rand_angle = []{ return rand()*1.0/RAND_MAX * 4. * M_PI - 2. * M_PI; };
    for (size_t i = 0; i < 99; i++) { 
      angles[i] = rand_angle(); 
    }

    for (size_t i = 0; i < 99; i++) {
        sin_results[i] = sin(angles[i]);
        fsine_results[i] = tscheb_fsine(angles[i], false);
        dsine_results[i] = tscheb_dsine(angles[i], false);
    }
    for (size_t i = 0; i < 99; i++) {
      angles_normalized[i] = pi_to_pi(angles[i]);
    }
    tscheb_dsines(angles_normalized, 99, dsine_vec_results);
    tscheb_dsines_unrolled(angles_normalized, 99, dsine_unrolled_results);

    for (size_t i = 0; i < 99; i++) {
      expect(that % fabs(fsine_results[i] - sin_results[i]) < 1e-6);
      expect(that % fabs(dsine_results[i] - sin_results[i]) < 1e-6);
      expect(that % fabs(dsine_vec_results[i] - sin_results[i]) < 1e-6);
      expect(that % fabs(dsine_unrolled_results[i] - sin_results[i]) < 1e-6);
    }

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&calc_tscheb_normalized_dsines)> bench("normalized_sine");
    bench.controls.NUM_RUNS = 50;
    bench.controls.REP = 5;
    bench.controls.CYCLES_REQUIRED = 0.0;
    bench.csv_path = "tscheb_angles.csv";
    // bench.csv_output = true;

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
#ifdef __AVX2__
    bench.add_function(&read_sine_vec, "read_sine_vec", 6*NR);
#endif

    // Run the benchmark: give the inputs of your function in the same order as they are defined. 
    for (size_t i = 5; i <= 16; i++) {
      NR = pow(2, i);
      std::cerr<< "Benchmarking N="<<NR<<" angles..."<<std::endl;

      bench.funcFlops[0] = NR*tp.sin; 
      bench.funcFlops[1] = NR*tscheb_dsine_flops(0, true);
      bench.funcFlops[2] = NR*tscheb_dsine_flops(0, true);
      bench.funcFlops[3] = tscheb_dsines_flops(alphas_d, NR, alphas_d);
      bench.funcFlops[4] = tscheb_dsines_unrolled_flops(alphas_d, NR, alphas_d);
      bench.funcFlops[5] = tscheb_dsines_avx_flops(alphas_d, NR, alphas_d);
      bench.funcFlops[6] = FlopCount::without_instr_mix(1);
      bench.funcFlops[7] = FlopCount::without_instr_mix(1);

      bench.funcBytes[0] = NR;
      bench.funcBytes[1] = NR;
      bench.funcBytes[2] = NR;
      bench.funcBytes[3] = NR;
      bench.funcBytes[4] = NR;
      bench.funcBytes[5] = NR;
      bench.funcBytes[6] = NR;
      bench.funcBytes[7] = NR;

      bench.run_name = std::to_string(NR); // Set name of run to identify it easier
      bench.run_benchmark(alphas_d);
    }

    bench.details();
    bench.write_csv_details();

    // Output is given when bench is destroyed (to change this behaviour, set bench.destructor_output=false). Optionally you can call bench.summary() to get it.

    // If you want the summary to be written to a file, set bench.fout with your preferred ostream. (Default std::cout)

    // Free memory
    free(angles);
    free(angles_normalized);
    destroy(alphas_d);

    return 0;
}
