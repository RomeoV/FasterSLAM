# FasterSLAM
> Fast Slam on steroids.

### Directory layout
Please consider this documented directory layout.
The structure is mostly preserved from the _yglee_ repository.
Note that the test directory has slightly changed.
More on tests below.
```
src
├── build/                          <-- run `cmake ..`, `make` and `ctest` here
├── CMakeLists.txt
├── core
│   ├── CMakeLists.txt              <-- add your source files here (no headers)
│   ├── example.cpp
│   ├── example.hpp
│   └── tests
│       ├── CMakeLists.txt
│       └── example_test.cpp        <-- testfiles ending in `_test.cpp` get registered automatically
├── external
│   ├── CMakeLists.txt
│   └── ut
│       └── ut.hpp
└── fastslam1
    ├── CMakeLists.txt              <-- add your source files here (no headers)
    ├── fastslam1.hpp               <-- header and source files with \
    ├── fastslam1.cpp               <-- the same name
    ├── main.cpp                    <-- main executable
    └── tests
        ├── CMakeLists.txt
        └── fastslam1_test.cpp      <-- testfiles ending in `_test.cpp`
```

### On testing
I have decided to use Kris Jusiak's lightweight testing framework [boost::ut](https://github.com/boost-experimental/ut) (Note that it is only a single header file and not part of the official boost library).
The library allows for a wide range of testing setups, including [Behaviour-Driven Development (BDD)](https://en.wikipedia.org/wiki/Behavior-driven_development) (see below for example) alongside "regular" tests.
I stronly encourage everyone to at least skim the excelent github page linked above.

The idea is that for every function/file `feature.hpp` we create a `feature_test.cpp` with a main and corresponding (Unit) Tests.
In there, we can use `boost::ut`'s testing functionality, as for example seen below in `example_test.cpp`.
Note that it is not neccessary to use the __BDD__ naming (_given_, _when_, _then_) but it makes for nicely structured tests.
```cpp
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
}
```

Note also that following the notion of [Extreme Programming](https://en.wikipedia.org/wiki/Extreme_programming) we might even follow the notion of defining the tests for a feature even before implementing the feature itself.
This way, the requirements are quite clear and can be well tested while developing/chaning other code.
