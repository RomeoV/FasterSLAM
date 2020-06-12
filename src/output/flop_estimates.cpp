#include <math.h>
#include <chrono>
#include <stdio.h>
#include <tsc_x86.h>
using namespace std;
using namespace std::chrono;

// timer cribbed from
// https://gist.github.com/gongzhitaao/7062087
class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return duration_cast<second_>
            (clock_::now() - beg_).count();
    }

private:
    typedef high_resolution_clock clock_;
    typedef duration<double, ratio<1>> second_;
    time_point<clock_> beg_;
};

int main(char* argv)
{
    double total;
    Timer tmr;
    init_tsc();
    double cycles;
    myInt64 start, end; 

#define randf() ((double) 3.0* rand()) / ((double) (RAND_MAX))
#define OP_TEST(name, expr)               \
    total = 0.0;                          \
    srand(42);                            \
    cycles = 0.0;                         \
    start = start_tsc() ;                 \
    for (int i = 0; i < 100000000; i++) { \
        double r1 = randf();              \
        double r2 = randf();              \
        total += expr;                    \
    }                                     \
    end = stop_tsc(start);                \
    double name = (double) end ;          \
                                          \
    printf(#name);                        \
    printf(" %.7f\n", (name - baseline)/322519040.0000000);

    // time the baseline code:
    //   for loop with no extra math op
    OP_TEST(baseline, 1.0)

    // time various floating point operations.
    //   subtracts off the baseline time to give
    //   a better approximation of the cost
    //   for just the specified operation

    OP_TEST(plus, r1 + r2)
    OP_TEST(mult, r1 * r2)
    OP_TEST(div, r1 / r2)
    OP_TEST(sqrt, sqrt(r1))
    OP_TEST(sin, sin(r1))
    OP_TEST(cos, cos(r1))
    OP_TEST(tan, tan(r1))
    OP_TEST(atan2, atan2(r1,r2))
    OP_TEST(exp, exp(r1))
    OP_TEST(greater, r1 > r2)
    OP_TEST(greaterequal, r1 >= r2)
    OP_TEST(random, rand())
    OP_TEST(abs, abs(r1))
    OP_TEST(floor, floor(r2))
    
    return 0;
}