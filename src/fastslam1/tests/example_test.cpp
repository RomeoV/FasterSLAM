#include "ut.hpp"
#include "example.hpp"
using namespace boost::ut;

int main() {
  "example test"_test = [] {
    expect(true);
  };
};
