#include "ut.hpp"
#include "nanobench.h"
#include "test_util.h"

#include <random>
#include <math.h>

#include "pi_to_pi.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

int main() {
    // Generate test data
    const size_t N = 1000;
    double angles[N];
    std::generate(angles, angles+N, []{return 2*M_PI*((rand()%1000)/100. - 5);});

    // Register functions to benchmark (with given name)
    auto test_functions = std::vector<NamedFunction<std::function<double (double)>>> {
        {"pi_to_pi", pi_to_pi},
        {"pi_to_pi_fmod", pi_to_pi_fmod},
    };

    // Compare two methods
    "functional equality"_test = [&] {
        auto is_close = [](auto lhs, auto rhs) -> bool {return lhs - rhs < 1e-14 and rhs - lhs < 1e-14;};
        for (size_t i = 0; i < N; i++) {
            expect(same_results(test_functions, is_close, angles[i]));
        }
    };

    // Benchmarking part
    static int i = 0;
    double sum = 0;
    auto bench = ankerl::nanobench::Bench();
    bench.minEpochIterations(100);
    bench.warmup(N);

    // Generic benchmarking test
    "benchmark methods"_test = [&] (auto f) {
        i = 0;
        sum = 0;
        bench.run(f.name, [&] {
            sum += f.function(angles[i%N]);
            i++;
        }).doNotOptimizeAway(sum);
    } | test_functions;
}
