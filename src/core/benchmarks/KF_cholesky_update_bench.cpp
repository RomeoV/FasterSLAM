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
#ifdef KF_YGLEE
    Vector2d x_init = {3.2403905331533212, -25.689432087069857};
    Matrix2d P_init = {0.20063369668655512, 0.018909593226709744, 0.018909593226709744, 0.011875705723671498};
#else
    double x_init[2] = {3.227460886446243, -25.613382543676146};
    double P_init[4] = {0.199855490073439, 0.019180472296076, 0.019180472296076, 0.011937739684843};
#endif
    std::copy(x_init, x_init+2, x);
    std::copy(P_init, P_init+4, P);
}

int main() {

    // Test: 
#ifdef KF_YGLEE
    cVector2d v = {-0.017001037783700212, -0.010645013219889199};
    cMatrix2d R = {0.010000000000000002, 0, 0, 0.00030461741909055634};
    cMatrix2d H = {0.073819203427568675, -0.99727164063023443, 0.038893076096335785, 0.0028789105989867826};    
#else
    double v[2] = {0.128762949830296, 0.019814250533567};
    double R[4] = {0.010000000000000, 0.0, 0.0, 0.000304617419787};
    double H[4] = {0.075431770172036, -0.997150965525639, 0.038902059641499, 0.002942835461779};
#endif
    Vector2d exact_x, x, x_v1/*, x_v2*/;
    Matrix2d exact_P, P, P_v1/*, P_v2*/;

    data_loader(exact_x, exact_P, v, R, H);
    KF_cholesky_update_base(exact_x, exact_P, v, R, H);
    
    data_loader(x_v1, P_v1, v, R, H);
    KF_cholesky_update_fused_ops(x_v1, P_v1, v, R, H);
    
    //data_loader(x_v2, P_v2, v, R, H);
    //KF_cholesky_update_v2(x_v2, P_v2, v, R, H);
     
    data_loader(x, P, v, R, H);
    KF_cholesky_update(x, P, v, R, H);

    // Check x
    double error = 0.0;
    for (int i = 0; i < 2; i++) {
        error = fabs(        x[i] - exact_x[i] ); expect(that % error < 1e-10) << i;
        error = fabs(     x_v1[i] - exact_x[i] ); expect(that % error < 1e-10) << i + 10;
        //error = fabs(     x_v2[i] - exact_x[i] ); expect(that % error < 1e-10) << i + 20;
    }
    // Check P 
    for (int i = 0; i < 4; i++) {
        error = fabs(        P[i] - exact_P[i] ); expect(that % error < 1e-10) << i;
        error = fabs(     P_v1[i] - exact_P[i] ); expect(that % error < 1e-10) << i + 10;
        //error = fabs(     P_v2[i] - exact_P[i] ); expect(that % error < 1e-10) << i + 20;
    }

    Benchmark<decltype(&KF_cholesky_update)> bench("KF_cholesky_update benchmark");
    double work = 120; // !!! Needs to be modified for certain variants !!!
    bench.data_loader = data_loader; // To guarantee same inputs
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&KF_cholesky_update_base, "base", work);
    bench.add_function(&KF_cholesky_update_fused_ops, "fused-ops", work);
    //bench.add_function(&KF_cholesky_update_v2, "v2", work);
    bench.add_function(&KF_cholesky_update, "active", work);

    bench.run_benchmark(x, P, v, R, H);

    return 0;
}

