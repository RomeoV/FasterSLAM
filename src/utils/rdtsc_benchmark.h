#pragma once

#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <iomanip>

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

    std::vector<comp_func> userFuncs;
    std::vector<std::string> funcNames;

    std::vector<int> funcFlops; //Work W
    std::vector<double> cycles; //Runtime T
    std::vector<double> speedups; //Speedup S
    std::vector<double> performances; //Performance P

    bool desctructor_output = true;

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

        if (!performances.empty()) {
            performances.clear();
            cycles.clear();
            speedups.clear();
        }

        for (int i = 0; i < numFuncs; i++)
        {
            double T = measure_cycles(userFuncs[i],std::forward<Args>(args)...);

            performances.push_back(funcFlops[i]/T);
            cycles.push_back(T);

            if (i== 0) {
                speedups.push_back(1.0);
            } else {
                speedups.push_back(cycles[0]/T);
            }
        }
    }

    void run_benchmark() {

        if (!performances.empty()) {
            performances.clear();
            cycles.clear();
            speedups.clear();
        }
        
        for (int i = 0; i < numFuncs; i++)
        {
            double T = measure_cycles(userFuncs[i]);

            performances.push_back(funcFlops[i]/T);
            cycles.push_back(T);

            if (i== 0) {
                speedups.push_back(1.0);
            } else {
                speedups.push_back(cycles[0]/T);
            }
        }
    }

    void summary(){
        assert(!performances.empty());

        fout<<name<<":"<<std::endl;
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
                                    <<separator<<prd(cycles[i], 4, cell_width) 
                                    <<separator<<prd(performances[i], 4, cell_width) 
                                    <<separator<<prd(speedups[i], 4, cell_width) 
                                    <<separator<<std::endl;
        }
    }

    ~Benchmark() {
        if (desctructor_output) {
            summary();
        }
    }
};

/* Global vars, used to keep track of student functions */


template<typename comp_func,
         typename ... Args>
double measure_cycles(comp_func f, Args&&... args) {
    double cycles = 0.;
    long num_runs = 100;
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
    long num_runs = 100;
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

