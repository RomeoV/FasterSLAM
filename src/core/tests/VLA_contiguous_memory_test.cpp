#include "ut.hpp"

using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
    given("A static sized array typedef") = [] {
        const size_t N_INNER = 33;
        const size_t N_OUTER = 100;
        typedef int FSA[N_INNER];  // fixed-size array

        "FSA is continuous"_test = [&] {
            FSA arr_1d;
            int* start = arr_1d;
            for (size_t i = 0; i < N_INNER; i++) {
                expect(that % &arr_1d[i] == start + i);
            }
        };

        when("I create another array of that") = [&] {
            FSA arr_2d[N_OUTER];
            then("I expect the memory to be continuous") = [&] {
                int* start = arr_2d[0];
                for (size_t outer = 0; outer < N_OUTER; outer++) {
                    for (size_t inner = 0; inner < N_INNER; inner++) {
                        expect(that % &arr_2d[outer][inner] == start + (N_INNER*outer + inner));
                    }
                }
            };
        };
    };
}