#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "data_associate_known.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

auto data_loader(cVector2d z[], const int* idz, const size_t idz_size, int* table, const int Nf_known, Vector2d zf[], int *idf, size_t *count_zf, Vector2d zn[], size_t *count_zn) {
    count_zf = 0;
    count_zn = 0;
    for (size_t i = 0; i < 35; i++) {
        table[i] = -1;
    }
}

int main() {

    // Test: 
    cVector2d z[2] = {{25.783699379426974, -1.4642259491817695}, 
                      {25.286434348683521, 0.14502450890426782}};
    const int idz[2] = {0, 21};

    Vector2d zf[35];     // will be modified in func
    Vector2d zn[35];     // will be modified in func
    int idf[35] = {};    // will be modified in func
    size_t count_zf = 0; // will be modified in func
    size_t count_zn = 0; // will be modified in func
    int table[35];       // will be modified in func
    int exact_table[35]; // will be modified in func
    
    data_loader(z, idz, 2, exact_table, 0, zf, idf, &count_zf, zn, &count_zn);
    // modifies table, zf, idf, zn and the count_zf, count_zh
    data_associate_known_base(z, idz, 2, exact_table, 0, zf, idf, &count_zf, zn, &count_zn);
    
    data_loader(z, idz, 2, table, 0, zf, idf, &count_zf, zn, &count_zn);
    // modifies table, zf, idf, zn and the count_zf, count_zh
    data_associate_known(z, idz, 2, table, 0, zf, idf, &count_zf, zn, &count_zn);
    
    // Check x
    for (int i = 0; i < 35; i++) {
        double error = fabs( table[i] - exact_table[i] );
        expect(that % error < 1e-12) << i;
    }

    Benchmark<decltype(&data_associate_known)> bench("data_associate_known benchmark");
    double work = 0; // TODO Count work
    bench.data_loader = data_loader; // To guarantee same inputs
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&data_associate_known_base, "base", work);
    bench.add_function(&data_associate_known, "active", work);

    bench.run_benchmark(z, idz, 2, table, 0, zf, idf, &count_zf, zn, &count_zn);

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

