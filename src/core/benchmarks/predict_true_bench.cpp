#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "predict_true.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

void data_loader(const double V,const double G,const double WB,
                const double dt, Vector3d xv) {
    xv[0] = 3.0;
    xv[1] = 1.0;
    xv[2] = 2.0;
}


int main() {

    double V = 1.0;
    double G = M_PI/3.0; //60°
    double WB = 0.1;
    double dt = 1.0;
    double xv[3];
    data_loader(V, G, WB, dt, xv); // Load xv

    double exact_xv[3];
    data_loader(V, G, WB, dt, exact_xv); // Load exact_xv

    for(int i = 0; i<3; i++){
        expect(that % fabs(xv[i]-exact_xv[i]) <= 1.0e-10) <<i ;
    }

    Benchmark<decltype(&predict_true)> bench("predict_true Benchmark");
    double work = 19;
    bench.data_loader = data_loader; // To guarantee same inputs
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&predict_true_base, "base", work);
    bench.add_function(&predict_true, "active", work);

    bench.run_benchmark(V, G, WB, dt, xv);

    G=0.0;
    bench.run_benchmark(V, G, WB, dt, xv);

    G=M_PI;
    bench.run_benchmark(V, G, WB, dt, xv);

    bench.destructor_output = false;

    bench.summary_long();

    return 0;
}