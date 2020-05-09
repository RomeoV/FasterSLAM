#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "KF_cholesky_update.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

auto data_loader(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H) {
    Vector2d x_init = {3.2403905331533212, -25.689432087069857};
    Matrix2d P_init = {0.20063369668655512, 0.018909593226709744, 0.018909593226709744, 0.011875705723671498};
    std::copy(x_init, x_init+2, x);
    std::copy(P_init, P_init+2, P);
}

int main() {

    // Test: 
    cVector2d v = {-0.017001037783700212, -0.010645013219889199};
    cMatrix2d R = {0.010000000000000002, 0, 0, 0.00030461741909055634};
    cMatrix2d H = {0.073819203427568675, -0.99727164063023443, 0.038893076096335785, 0.0028789105989867826};    
    Vector2d exact_x, x, x_v1, x_v2;
    Matrix2d exact_P, P, P_v1, P_v2;
    
    data_loader(exact_x, exact_P, v, R, H);
    KF_cholesky_update_base(exact_x, exact_P, v, R, H);
    
    data_loader(x_v1, P_v1, v, R, H);
    KF_cholesky_update_v1(x_v1, P_v1, v, R, H);
    
    data_loader(x_v2, P_v2, v, R, H);
    KF_cholesky_update_v2(x_v2, P_v2, v, R, H);
    
    data_loader(x, P, v, R, H);
    KF_cholesky_update(x, P, v, R, H);

    // Check x
    double error = 0.0;
    for (int i = 0; i < 2; i++) {
        error = fabs(    x[i] - exact_x[i] ); expect(that % error < 1e-12) << i;
        error = fabs( x_v1[i] - exact_x[i] ); expect(that % error < 1e-12) << i + 10;
        error = fabs( x_v2[i] - exact_x[i] ); expect(that % error < 1e-12) << i + 20;
    }
    // Check P 
    for (int i = 0; i < 4; i++) {
        error = fabs(    P[i] - exact_P[i] ); expect(that % error < 1e-12) << i;
        error = fabs( P_v1[i] - exact_P[i] ); expect(that % error < 1e-12) << i + 10;
        error = fabs( P_v2[i] - exact_P[i] ); expect(that % error < 1e-12) << i + 20;
    }

    Benchmark<decltype(&KF_cholesky_update)> bench("KF_cholesky_update benchmark");
    double work = 120; // !!! Needs to be modified for certain variants !!!
    bench.data_loader = data_loader; // To guarantee same inputs
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&KF_cholesky_update_base, "base", work);
    bench.add_function(&KF_cholesky_update_v1, "v1", work);
    bench.add_function(&KF_cholesky_update_v2, "v2", work);
    bench.add_function(&KF_cholesky_update, "active", work);

    bench.run_benchmark(x, P, v, R, H);

    /*
    // Alternative (much slower here, but nicer to look at. Generally useful if you want to average over a few inputs). Yields averages over all runs.
    
    Benchmark<decltype(&pi_to_pi)> bench("pi_to_pi Benchmark");

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&pi_to_pi, "pi_to_pi", 6);
    bench.add_function(&pi_to_pi_fmod, "pi_to_pi_fmod", 6);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    for (int i = 0; i<N; i++) {
        // You could set the data_loader function here to generate new input. bench.data_loader =&my_load_func_i...
        bench.run_benchmark(angles[i]);
    }
    */

    return 0;
}

