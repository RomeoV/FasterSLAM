#include "get_observations.h"  // import file to test
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

    "get_observations"_test = [] {
        given("I have some arguments") = [] {
            const double rmax = 30; 
            const double x[3] = {0.674090417131751, -0.030904471130924017, -0.0073589032333721818};
            double z[4] = { };
            size_t lm_cols = 35, nidf = 35;

            FILE* fp = fopen("../../../core/tests/lm.txt", "r"); 

            double *lm = (double*) malloc( 2*lm_cols * sizeof(double) ); 
            for (size_t i = 0; i < lm_cols; i++) {
                fscanf(fp, "%lf\t%lf\n", &lm[0*lm_cols+i], &lm[1*lm_cols+i]);
            }

            int *idf = (int*) malloc( nidf * sizeof(int) );
            for (int i = 0; i < nidf; i++) {
                idf[i] = i;
            }

            when("I call get_observations()") = [&] {

                get_observations(x, rmax, lm, lm_cols, &idf, &nidf, z);

                then("This is equal with the actual result") = [&] {
                    const double actual_z[4] = {25.77444969441645, -1.4733774573801159, 25.276107769232738, 0.13836675004531551};
                    for (int i = 0; i < 4; i++) {    
                        expect(fabs(z[i] - actual_z[i]) < 1e-06) << std::setprecision(12) << z[i] << " != " << actual_z[i];
                    }
                };
            };

            fclose( fp );
            free( lm );
            free( idf );
        };
    };
}