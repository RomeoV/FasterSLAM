#include "ut.hpp"
#include "nanobench.h"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "pi_to_pi.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

int main() {
    // Generate test data
    const size_t N = 1000;
    double angles[N];
    std::generate(angles, angles+N, []{return 2*M_PI*((rand()%1000)/100. - 5);});

    // Compare two methods
    "functional equality"_test = [&] {
        auto is_close = [](auto lhs, auto rhs) -> bool {return lhs - rhs < 1e-14 or rhs - lhs < 1e-14;};
        for (size_t i = 0; i < N; i++) {
            double method1 = pi_to_pi(angles[i]);
            double method2 = pi_to_pi_fmod(angles[i]);
            expect(is_close(method1, method2));
        }
    };

    // Benchmarking part
    static int i = 0;
    double sum = 0;
    auto bench = ankerl::nanobench::Bench();
    bench.minEpochIterations(100);
    bench.warmup(N);

    // Maybe this could be a bit clearer
    // Register functions to benchmark (with given name)
    using TestFunctionSignature = std::function<double (double)>;
    using NamedTestFunctions = std::vector<std::pair<std::string, TestFunctionSignature>>;
    NamedTestFunctions test_functions = {{"pi_to_pi", pi_to_pi},
                                         {"pi_to_pi_fmod", pi_to_pi_fmod}};

    // Generic benchmarking test
    "benchmark methods"_test = [&] (auto f) {
        i = 0;
        bench.run(f.first, [&] {
            sum += f.second(angles[i%N]);
            i++;
        }).doNotOptimizeAway(sum);
    } | test_functions;
}
