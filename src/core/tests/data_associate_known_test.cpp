#include "data_associate_known.h"  // import file to test
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

    "data_associate_known"_test = [] {
        given("I have some arguments") = [] {

            cVector2d z[2] = {{25.783699379426974, -1.4642259491817695}, 
                               {25.286434348683521, 0.14502450890426782}};
            const int idz[2] = {0, 21};
            Vector2d zf[35]; // will be modified in func
            Vector2d zn[35]; // will be modified in func
            int idf[35] = {};   // will be modified in func
            size_t count_zf = 0;
            size_t count_zn = 0;

            int table[35]; // initilize here, will be modified in func
            for (size_t i = 0; i < 35; i++) {
                table[i] = -1;
            }

            when("I call data_associate_known()") = [&] {
                int table_target[35];
                for (size_t  i = 0; i < 35; i++) {
                    table_target[i] = -1;
                }
                table_target[0] = 0;
                table_target[21] = 1;
                
                const Vector2d zn_target[2] = {25.783699379426974, -1.4642259491817695, 
                                             25.286434348683521, 0.14502450890426782};

                // modifies table, zf, idf, zn and the count_zf, count_zh
                data_associate_known(z, idz, 2, table, 0, zf, idf, &count_zf, zn, &count_zn);

                then("This is equal with the actual result [table]") = [&] {
                    for (size_t i = 0; i < 35; i++) {
                        expect(abs(table[i] - table_target[i]) < 1e-10) << std::setprecision(12) << table[i] << " != " << table_target[i];
                    }
                };
                // zn needs to be allocated inside data_associate_known()
                then("This is equal with the actual result [zn]") = [&] {
                    for (size_t i = 0; i < 2; i++) {
                        double zni0 = zn[i][0];
                        double zni1 = zn[i][1];
                        double zn_target_i0 = zn[i][0];
                        double zn_target_i1 = zn[i][1];
                        expect(fabs(zni0 - zn_target_i0) < 1e-10) << std::setprecision(12) << zni0 << " != " << zn_target_i0;
                        expect(fabs(zni1 - zn_target_i1) < 1e-10) << std::setprecision(12) << zni1 << " != " << zn_target_i1;
                    }
                };
                // zf, idf still NULL in this case
            };
        };
    };
}
