#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <iostream>
#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>
#include <immintrin.h>

#include "trigonometry.h"
#include "linalg.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`



int main() {

    // Compare two methods
    /*"functional equality"_test = [&] {
        auto is_close = [&](double lhs, double rhs) -> bool {return std::abs(lhs - rhs) < 1e-4;};
        expect(is_close(Hf_base[i][j], Hf[i][j])) <<i << "Hf" <<j;
    };*/

    auto ymm0 = _mm256_set_pd(3,2,1,0);

    // no functional equality test
    Benchmark<decltype(&atan2)> bench("Atan2 Benchmark");
    bench.add_function(&atan2_approximation1, "atan2_approx1", 0.0);
    bench.run_benchmark(1.0, 0.0);
#ifdef __AVX2__
    bench.add_function(&atan2_approximation2, "atan2_approx2", 0.0);
    bench.run_benchmark(1.0, 0.0);
#endif

    bench.details();

    return 0;
}






 
