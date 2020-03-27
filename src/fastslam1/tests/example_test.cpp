#include "ut.hpp"
#include "example.hpp"
using namespace boost::ut;

int main() {
  "example test"_test = [] {
    expect(true);
    expect(sum(1,2,3) == 6) << "Check the sum method for some argument";
  };
};
