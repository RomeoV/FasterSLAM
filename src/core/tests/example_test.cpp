#include "example.hpp"  // import file to test

#include <vector>  // used for test input
#include "ut.hpp"  // import functionality and namespaces from single header file
using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

int main() {
  "example test"_test = [] {
    expect(true) << "with more information!";
  };

  "vector add"_test = [] {
    given("I have two arrays") = [] {
      const size_t N = 5;
      double lhs[N] = {1,2,3,4,5};
      double rhs[N] = {1,2,3,4,5};

      when("I add them") = [&] {
        double res[N];
        add_two_arrays(&lhs[0], &rhs[0], &res[0], N);

        then("I get the elementwise sum") = [=] {
          "elementwise equal"_test = [res](size_t i) {
            expect(res[i-1] == 2*i);
          } | std::vector{1,2,3,4,5};
        };
      };
    };
  };
};
