#include "rdtsc_benchmark.h"
#include <iostream>

#define NR 32

void rands(double * m, size_t row, size_t col)
{
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (size_t i = 0; i < row*col; ++i)  
        m[i] = dist(gen);
}

void build(double **a, int m, int n)
{
    *a = static_cast<double *>(aligned_alloc(32, m * n * sizeof(double))); //TODO: align this to 32
    rands(*a, m, n);
}

void destroy(double * m)
{
    free(m);
}

void kernel_base(double* A, double *x, double *y) {
  for (int i = 0; i < NR; ++i) {
    for (int j = 0; j < NR; ++j) {
      y[i] += (j + 2) * A[i*NR + j] * x[j];
    }
  }
  
  for (int i = 0; i < NR; ++i) {
    if (y[i] >= 0) {
      y[i] += 1;
    }
  }
}


void kernel_fast(double* A, double *x, double *y) {
  for (int i = 0; i < 16; ++i) { //Simple showcase: Just want lower runtime
    for (int j = 0; j < NR; ++j) {
      y[i] += (j + 2) * A[i*NR + j] * x[j];
    }
  }
  
  for (int i = 0; i < NR; ++i) {
    if (y[i] >= 0) {
      y[i] += 1;
    }
  }
}


int main() {
    // Initialize Input
    double *A, *x, *y;
    int n = NR;

    build(&A, n, n);
    build(&x, 1, n);
    build(&y, 1, n);

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&kernel_base)> bench("Example Benchmark");

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&kernel_base, "Kernel base", 3*NR*NR+NR);
    bench.add_function(&kernel_fast, "Kernel fast", 3*NR*NR+NR);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(A,x,y);

    // Output is given when bench is destroyed (to change this behaviour, set bench.destructor_output=false). Optionally you can call bench.summary() to get it.

    //If you want the summary to be written to a file, set bench.fout with your preferred ostream. (Default std::cout)

    // Free memory
    destroy(A);
    destroy(x);
    destroy(y);

    return 0;
}
