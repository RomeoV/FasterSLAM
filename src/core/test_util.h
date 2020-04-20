#pragma once

#include <functional>
#include <algorithm>
#include <numeric>
#include <vector>
#include <variant>
#include <string>

/** Function wrapper with name attached, mostly for microbenchmark output */
template <typename F>
struct NamedFunction {
    std::string name;
    F function;
};

struct IsEqual {
    bool is_equal;
    std::string error_msg;
};

/** Compares vector of named functions for result equality
 *
 *  @param test_functions Vector of NamedFunction
 *  @param is_close Function that should take two result types of the test_function and return true if they are almost equal
 *  @param args The arguments to be evaluated with the test functions 
 */
template <typename Func, typename Comp, typename ... Args>
auto same_results(std::vector<NamedFunction<Func>> test_functions, Comp is_close, Args... args) -> IsEqual {
    using result_type = typename Func::result_type;
    auto results = std::vector<result_type>(test_functions.size());
    std::transform(test_functions.begin(), test_functions.end(),  // input
                   results.begin(),  // output
                   [&](auto named_function){
                     return named_function.function(args...);
                   }
    );
    auto adjacent_equal = std::vector<bool>(test_functions.size());
    std::adjacent_difference(results.begin(), results.end(),  // input
                             adjacent_equal.begin(),  // output
                             [is_close](auto lhs, auto rhs) {
                                return is_close(lhs, rhs);
                             }
    );
    bool all_equal =  std::reduce(std::next(adjacent_equal.begin()), adjacent_equal.end(),
                                  true,
                                  std::logical_and{}
    );
    if (all_equal) {
        return {true, ""};
    }
    else {
        auto not_eq_it = std::find(std::next(adjacent_equal.begin()), adjacent_equal.end(), false);
        size_t dist = std::distance(adjacent_equal.begin(), not_eq_it);
        std::string lhs_func_name = test_functions.at(dist-1).name;
        std::string rhs_func_name = test_functions.at(dist  ).name;
        return {false, lhs_func_name + " != " + rhs_func_name};
    }
}
