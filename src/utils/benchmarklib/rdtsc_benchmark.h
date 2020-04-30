#pragma once

#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <iomanip>
#include <limits>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "string_format.h"
#include "tsc_x86.h"
#include <assert.h> 

//Macros (overwrite them in the benchmark if necessary)
#define CYCLES_REQUIRED 1e8
#define REP 50
#define NUM_RUNS 100

//headers
template<typename comp_func,
         typename ... Args>
double measure_cycles(comp_func f, Args&&... args);

template<typename comp_func>
double measure_cycles(comp_func f);

//Main struct
template <typename comp_func>
struct Benchmark {
    std::string name = "";

    // Optional function that takes the same input as comp_func, and resets the data.
    // Malloc should be done outside this function, here we only set values (fill arrays).
    // We do this for comparability of benchmarks where funcs overwrite the input.
    comp_func data_loader = nullptr;

    // comp_func work_compute = nullptr; //returns a vector that overwrites funcFlops, WIP

    std::vector<comp_func> userFuncs;
    std::vector<std::string> funcNames;

    std::vector<int> funcFlops; //Work W
    std::vector<double> cycles; // Sum of runtimes over all runs T
    std::vector<double> min_cycles; // Minimum t over all runs
    std::vector<double> max_cycles; // Maximum t over all runs
    std::vector<double> speedups; //Speedup S
    std::vector<double> performances; //Performance P

    int num_runs = 0;

    bool destructor_output = true;

    std::ostream & fout = std::cout;
    int numFuncs = 0;
    Benchmark(std::string title) {
        name = title;
    }
    void add_function(comp_func f, std::string name, int flop){
        userFuncs.push_back(f);
        funcNames.emplace_back(name);
        funcFlops.push_back(flop);
        numFuncs++;
    }

    template <typename ... Args>
    void run_benchmark(Args&&... args) {
        
        
        num_runs++;
        if (performances.empty()) {
            performances.resize(numFuncs, 0.0);
            cycles.resize(numFuncs, 0.0);
            min_cycles.resize(numFuncs, std::numeric_limits<float>::max());
            max_cycles.resize(numFuncs, -std::numeric_limits<float>::max());
            speedups.resize(numFuncs, 0.0);
            speedups[0] = 1.0 ; //Reference function at index 0
        }

        
        for (int i = 0; i < numFuncs; i++)
        {
            if (data_loader) {
                data_loader(std::forward<Args>(args)...);
            }
            
            double t = measure_cycles(userFuncs[i],std::forward<Args>(args)...);

            double T = cycles[i] + t;

            min_cycles[i] = std::min(min_cycles[i],t);
            max_cycles[i] = std::max(max_cycles[i],t);

            cycles[i] = T;
            
            performances[i] = funcFlops[i] * double(num_runs) / T ;

            if (i > 0) {
                speedups[i] = cycles[0]/T;
            }
        }
    }

    void run_benchmark() {
        num_runs++;
        if (performances.empty()) {
            performances.resize(numFuncs, 0.0);
            cycles.resize(numFuncs, 0.0);
            min_cycles.resize(numFuncs, std::numeric_limits<float>::max());
            max_cycles.resize(numFuncs, -std::numeric_limits<float>::max());
            speedups.resize(numFuncs, 0.0);
            speedups[0] = 1.0 ; //Reference function at index 0
        }

        
        for (int i = 0; i < numFuncs; i++)
        {
            if (data_loader) {
                data_loader();
            }
            
            double t = measure_cycles(userFuncs[i]);

            double T = cycles[i] + t;

            min_cycles[i] = std::min(min_cycles[i],t);
            max_cycles[i] = std::max(max_cycles[i],t);

            cycles[i] = T;
            
            performances[i] = funcFlops[i] * double(num_runs) / T ;

            if (i > 0) {
                speedups[i] = cycles[0]/T;
            }
        }
    }

    void summary(){
        assert(!performances.empty());

        fout<<name<<" ("<<num_runs<<" runs, avg):"<<std::endl;
        const int cell_width = 17;
        const char* separator = " | ";
        const char* underline = "\033[4m";
        const char* no_underline = "\033[0m";
        fout<<underline<<separator<<right("i",2)
                            <<separator<<left("Name",30)
                            <<separator<<right("Work [flops]",cell_width) 
                            <<separator<<right("Time [cyc]",cell_width) 
                            <<separator<<right("Perf. [flops/cyc]",cell_width) 
                            <<separator<<right("Speedup [-]",cell_width) 
                            <<separator<<std::endl;
        
        for (int i = 0; i<numFuncs;i++){
            fout<<no_underline<<separator<<prd(i, 0,2)
                                    <<separator<<left(funcNames[i],30)
                                    <<separator<<prd(funcFlops[i], 0, cell_width)
                                    <<separator<<prd(cycles[i] / num_runs, 4, cell_width) 
                                    <<separator<<prd(performances[i], 4, cell_width) 
                                    <<separator<<prd(speedups[i], 4, cell_width) 
                                    <<separator<<std::endl;
        }
    }

    void summary_long(){
        assert(!performances.empty());

        fout<<name<<" ("<<num_runs<<" runs, avg):"<<std::endl;
        const int cell_width = 17;
        const char* separator = " | ";
        const char* underline = "\033[4m";
        const char* no_underline = "\033[0m";
        fout<<underline<<separator<<right("i",2)
                            <<separator<<left("Name",30)
                            <<separator<<right("Work [flops]",cell_width) 
                            <<separator<<right("Time [cyc]",cell_width) 
                            <<separator<<right("min. t [cyc]",cell_width)
                            <<separator<<right("max. t [cyc]",cell_width)
                            <<separator<<right("Perf. [flops/cyc]",cell_width) 
                            <<separator<<right("Speedup [-]",cell_width) 
                            <<separator<<std::endl;
        
        for (int i = 0; i<numFuncs;i++){
            fout<<no_underline<<separator<<prd(i, 0,2)
                                    <<separator<<left(funcNames[i],30)
                                    <<separator<<prd(funcFlops[i], 0, cell_width)
                                    <<separator<<prd(cycles[i] / num_runs, 4, cell_width)
                                    <<separator<<prd(min_cycles[i], 4, cell_width) 
                                    <<separator<<prd(max_cycles[i], 4, cell_width)  
                                    <<separator<<prd(performances[i], 4, cell_width) 
                                    <<separator<<prd(speedups[i], 4, cell_width) 
                                    <<separator<<std::endl;
        }
    }

    ~Benchmark() {
        if (destructor_output) {
            summary();
        }
    }
};

/* Global vars, used to keep track of student functions */


template<typename comp_func,
         typename ... Args>
double measure_cycles(comp_func f, Args&&... args) {
    double cycles = 0.;
    long num_runs = NUM_RUNS;
    double multiplier = 1;
    myInt64 start, end;

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            std::forward<comp_func>(f)(std::forward<Args>(args)...);          
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);
        
    } while (multiplier > 2);

    std::list<double> cyclesList;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            std::forward<comp_func>(f)(std::forward<Args>(args)...);            
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cyclesList.push_back(cycles);
    }
    total_cycles /= REP;

    cyclesList.sort();
    cycles = total_cycles;
    return  cycles;
}

template<typename comp_func>
double measure_cycles(comp_func f) {
    double cycles = 0.;
    long num_runs = NUM_RUNS;
    double multiplier = 1;
    myInt64 start, end;

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            f();          
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);
        
    } while (multiplier > 2);

    std::list<double> cyclesList;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            f();            
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cyclesList.push_back(cycles);
    }
    total_cycles /= REP;

    cyclesList.sort();
    cycles = total_cycles;
    return  cycles;
}

