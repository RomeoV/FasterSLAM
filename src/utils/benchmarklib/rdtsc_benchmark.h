#pragma once

#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "string_format.h"
#include "tsc_x86.h"
#include "flop_count.h"
#include <assert.h> 

//Macros 
#define CYCLES_REQUIRED_ 1e8 // Cache warmup
#define REP_ 50 // repetitions
#define NUM_RUNS_ 100 //number of runs within each repetition (avg)
#define WARMUP

typedef struct BenchControls {
    long CYCLES_REQUIRED = CYCLES_REQUIRED_;
    long REP = REP_;
    long NUM_RUNS = NUM_RUNS_;
} BenchControls;

//headers
template<typename comp_func,
         typename ... Args>
double measure_cycles(comp_func f, comp_func data_loader, BenchControls controls, Args&&... args);

//headers
template<typename comp_func,
         typename ... Args>
double measure_cycles(comp_func f, BenchControls controls, Args&&... args);

template<typename comp_func>
double measure_cycles(comp_func f, BenchControls controls);

template<typename comp_func>
double measure_cycles(comp_func f);

template<typename comp_func,
         typename ... Args>
double measure_cycles(comp_func f, Args&&... args);

//Main struct
template <typename comp_func>
struct Benchmark {
    std::string name = "";

    BenchControls controls; // To control repeats in benchmark!

    // Optional function that takes the same input as comp_func, and resets the data.
    // Malloc should be done outside this function, here we only set values (fill arrays).
    // We do this for comparability of benchmarks where funcs overwrite the input.
    comp_func data_loader = nullptr;

    // comp_func work_compute = nullptr; //returns a vector that overwrites funcFlops, WIP

    std::vector<comp_func> userFuncs;
    std::vector<std::string> funcNames;

    int numFuncs = 0;
    std::vector<FlopCount> funcFlops; //Work W
    std::vector<double> funcBytes; // Memory Moved Q [Bytes!]
    std::vector<double> flops_sum; //Sum of Work W for each function over all runs
    std::vector<double> bytes_sum; //Sum of Memory Moved Q for each function over all runs
    std::vector<double> cycles; // Sum of runtimes over all runs T
    std::vector<double> min_cycles; // Minimum t over all runs
    std::vector<double> max_cycles; // Maximum t over all runs
    std::vector<double> speedups; //Speedup S
    std::vector<double> performances; //Performance P
    
    std::vector<std::string> run_names;

    std::vector<std::vector<double>> cycles_capture;
    std::vector<std::vector<FlopCount>> flops_capture;
    std::vector<std::vector<double>> bytes_capture;
    
    int num_runs = 0;
    std::string run_name = std::to_string(num_runs);
    bool destructor_output = true;

    std::ostream & fout = std::cout;

    bool csv_output = true;
    std::string csv_path = "benchmark.csv";
    
    
    Benchmark(std::string title) {
        name = title;
    }

    void add_function(comp_func f, std::string name, FlopCount flop){
        double _nan = nan("1");
        add_function(f, name, flop, _nan);
    }

    void add_function(comp_func f, std::string name, double flop){
        add_function(f, name, FlopCount::without_instr_mix(flop));
    }

    void add_function(comp_func f, std::string name, FlopCount flop, double bytes){
        userFuncs.push_back(f);
        funcNames.emplace_back(name);
        funcFlops.push_back(flop);
        funcBytes.push_back(bytes);
        numFuncs++;
    }

    template <typename ... Args>
    void run_benchmark(Args&&... args) {
        
        run_names.push_back(run_name);
        
        if (performances.empty()) {
            performances.resize(numFuncs, 0.0);
            cycles.resize(numFuncs, 0.0);
            flops_sum.resize(numFuncs,0.0);
            bytes_sum.resize(numFuncs,0.0);
            cycles_capture.resize(numFuncs,{});
            flops_capture.resize(numFuncs,{});
            bytes_capture.resize(numFuncs,{});
            speedups.resize(numFuncs, 0.0);
            speedups[0] = 1.0 ; //Reference function at index 0
        }

        
        for (int i = 0; i < numFuncs; i++)
        {
            double t = measure_cycles(userFuncs[i], data_loader, controls,
                                                std::forward<Args>(args)...);

            double T = cycles[i] + t;

            cycles_capture[i].push_back(t);
            cycles[i] = T;

            flops_capture[i].push_back(funcFlops[i]);
            flops_sum[i] += funcFlops[i].flop_sum;

            bytes_capture[i].push_back(funcBytes[i]);
            bytes_sum[i] += funcBytes[i];
            
            performances[i] = flops_sum[i] / T ;

            if (i > 0) {
                speedups[i] = cycles[0]/T;
            }
        }
        num_runs++;
        run_name = std::to_string(num_runs);
        
    }

    void run_benchmark() {
        
        run_names.push_back(run_name);
        
        if (performances.empty()) {
            performances.resize(numFuncs, 0.0);
            cycles.resize(numFuncs, 0.0);
            flops_sum.resize(numFuncs,0.0);
            bytes_sum.resize(numFuncs,0.0);
            cycles_capture.resize(numFuncs,{});
            flops_capture.resize(numFuncs,{});
            bytes_capture.resize(numFuncs,{});
            speedups.resize(numFuncs, 0.0);
            speedups[0] = 1.0 ; //Reference function at index 0
        }

        
        for (int i = 0; i < numFuncs; i++)
        {
            if (data_loader) {
                data_loader();
            }

            double t = measure_cycles(userFuncs[i], controls);

            double T = cycles[i] + t;

            cycles_capture[i].push_back(t);
            cycles[i] = T;

            flops_capture[i].push_back(funcFlops[i]);
            flops_sum[i] += funcFlops[i].flop_sum;

            bytes_capture[i].push_back(funcBytes[i]);
            bytes_sum[i] += funcBytes[i];

            performances[i] = flops_sum[i] / T ;

            if (i > 0) {
                speedups[i] = cycles[0]/T;
            }
        }
        num_runs++;
        run_name = std::to_string(num_runs);
        
    }


    // Summarized output
    void summary(){
        assert(!performances.empty());
        

        fout<<name<<" ("<<num_runs<<" runs, avg):"<<std::endl;
        const int cell_width = 17;
        const char* separator = " | ";
        const char* underline = "\033[4m";
        const char* no_underline = "\033[0m";
        fout<<underline<<separator<<right("i",2)
                            <<separator<<left("Name",40)
                            <<separator<<right("Work [flops]",cell_width) 
                            <<separator<<right("Memory [bytes]",cell_width) 
                            <<separator<<right("Time [cyc]",cell_width) 
                            <<separator<<right("Perf. [flops/cyc]",cell_width) 
                            <<separator<<right("Speedup [-]",cell_width) 
                            <<separator<<std::endl;
        
        for (int i = 0; i<numFuncs;i++){
            fout<<no_underline<<separator<<prd(i, 0,2)
                                    <<separator<<left(funcNames[i],40)
                                    <<separator<<prd(flops_sum[i] / num_runs, 0, cell_width)
                                    <<separator<<prd(bytes_sum[i] / num_runs, 0, cell_width)
                                    <<separator<<prd(cycles[i] / num_runs, 4, cell_width) 
                                    <<separator<<prd(performances[i], 4, cell_width) 
                                    <<separator<<prd(speedups[i], 4, cell_width) 
                                    <<separator<<std::endl;
        }
    }

    // Summarized output, but with some extra columns
    void summary_long(){
        assert(!performances.empty());

        fout<<name<<" ("<<num_runs<<" runs, avg):"<<std::endl;
        const int cell_width = 17;
        const char* separator = " | ";
        const char* underline = "\033[4m";
        const char* no_underline = "\033[0m";
        fout<<underline<<separator<<right("i",2)
                            <<separator<<left("Name",40)
                            <<separator<<right("Work [flops]",cell_width) 
                            <<separator<<right("Memory [bytes]",cell_width) 
                            <<separator<<right("Time [cyc]",cell_width) 
                            <<separator<<right("min. t [cyc]",cell_width)
                            <<separator<<right("max. t [cyc]",cell_width)
                            <<separator<<right("Perf. [flops/cyc]",cell_width) 
                            <<separator<<right("Speedup [-]",cell_width) 
                            <<separator<<std::endl;
        
        for (int i = 0; i<numFuncs;i++){
            fout<<no_underline<<separator<<prd(i, 0,2)
                                    <<separator<<left(funcNames[i],40)
                                    <<separator<<prd(flops_sum[i] / num_runs, 0, cell_width)
                                    <<separator<<prd(bytes_sum[i] / num_runs, 0, cell_width)
                                    <<separator<<prd(cycles[i] / num_runs, 4, cell_width)
                                    <<separator<<prd(*std::min_element(cycles_capture[i].begin(),
                                                     cycles_capture[i].end()), 4, cell_width) 
                                    <<separator<<prd(*std::max_element(cycles_capture[i].begin(),
                                                     cycles_capture[i].end()), 4, cell_width)  
                                    <<separator<<prd(performances[i], 4, cell_width) 
                                    <<separator<<prd(speedups[i], 4, cell_width) 
                                    <<separator<<std::endl;
        }
    }

    // Outputs per func and run
    void details() {
        
        assert(!performances.empty());
        
        int run_name_size = 3;
        for (int j = 0; j<num_runs; j++) {
            run_name_size = std::max((int)run_names[j].length(), run_name_size);
        }
        fout<<name<<" ("<<num_runs<<" runs, avg):"<<std::endl;
        const int cell_width = 17;
        const char* separator = " | ";
        const char* underline = "\033[4m";
        const char* no_underline = "\033[0m";
        fout<<underline<<separator<<right("i",2)
                            <<separator<<left("Name",40)
                            <<separator<<right("Run",run_name_size)
                            <<separator<<right("Work [flops]",cell_width) 
                            <<separator<<right("Memory [bytes]",cell_width) 
                            <<separator<<right("Time [cyc]",cell_width) 
                            <<separator<<right("Perf. [flops/cyc]",cell_width) 
                            <<separator<<right("Gap [cyc/flop]",cell_width)  //Cheat
                            <<separator<<right("Speedup [-]",cell_width) 
                            <<separator<<std::endl;

        for (int i = 0; i<numFuncs;i++){
            for (int j = 0; j<num_runs; j++) {
                double cyc = cycles_capture[i][j];
                double flops = flops_capture[i][j].flop_sum;
                double bytes = bytes_capture[i][j];
                fout<<no_underline<<separator<<((j==0) ? prd(i, 0,2) : "  ")
                        <<separator<<left((j==0) ? funcNames[i] : "  ",40)
                        <<separator<<right(run_names[j],run_name_size)
                        <<separator<<prd(flops, 0, cell_width)
                        <<separator<<prd(bytes, 0, cell_width)
                        <<separator<<prd(cyc, 4, cell_width) 
                        <<separator<<prd(flops/cyc, 4, cell_width) 
                        <<separator<<prd(cyc/flops, 4, cell_width) 
                        <<separator<<prd(cycles_capture[0][j]/cyc, 4, cell_width) 
                        <<separator<<std::endl;
            }
        }
    }

    void write_csv() {
        std::ofstream fstream;
        fstream.open(csv_path, std::ios::out | std::ios::app);
        const char* separator = ";";
        const int cell_width = 0;
        for (int i = 0; i<numFuncs;i++){
            fstream<<name
                    <<separator<<left(funcNames[i], cell_width)
                    <<separator<<right("",cell_width)
                    <<separator<<prd(flops_sum[i] / num_runs, 0,  cell_width)
                    <<separator<<prd(bytes_sum[i] / num_runs, 0,  cell_width)
                    <<separator<<prd(cycles[i] / num_runs, 4,  cell_width) 
                    <<separator<<prd(performances[i], 4,  cell_width) 
                    <<separator<<prd(speedups[i], 4,  cell_width) 
                    <<std::endl;
        }
        fstream.close();
    }

    void write_csv_details() {
        std::ofstream fstream;
        fstream.open(csv_path, std::ios::out | std::ios::app);
        const char* separator = ";";
        const int cell_width = 0;

        for (int i = 0; i<numFuncs;i++){
            for (int j = 0; j<num_runs; j++) {
                double cyc = cycles_capture[i][j];
                double flops = flops_capture[i][j].flop_sum;
                double bytes = bytes_capture[i][j];
                fstream<<name<<separator<<left(funcNames[i], cell_width)
                        <<separator<<right(run_names[j],cell_width)
                        <<separator<<prd(flops, 0, cell_width)
                        <<separator<<prd(bytes, 0, cell_width)
                        <<separator<<prd(cyc, 4, cell_width) 
                        <<separator<<prd(flops/cyc, 4, cell_width) 
                        <<separator<<prd(cycles_capture[0][j]/cyc, 4, cell_width) 
                        <<std::endl;
            }
        }
        fstream.close();
    }

    ~Benchmark() {
        if (destructor_output) {
            summary();
        }
        if (csv_output) {
            write_csv();
        }
    }

    
};

/* Global vars, used to keep track of student functions */


template<typename comp_func,
         typename ... Args>
double measure_cycles(comp_func f, comp_func data_loader, BenchControls controls, Args&&... args) {
    double cycles = 0.;
    long num_runs = controls.NUM_RUNS;
    double multiplier = 1;
    myInt64 start, end;

    if (data_loader) {
        std::forward<comp_func>(data_loader)(std::forward<Args>(args)...);
    }
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
#ifdef WARMUP
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            std::forward<comp_func>(f)(std::forward<Args>(args)...);          
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (controls.CYCLES_REQUIRED) / (cycles);
        
    } while (multiplier > 2);
#endif

    std::list<double> cyclesList;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    
    double total_cycles = 0;
    for (size_t j = 0; j < controls.REP; j++) {

        if (data_loader) {
            std::forward<comp_func>(data_loader)(std::forward<Args>(args)...);
        }

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            std::forward<comp_func>(f)(std::forward<Args>(args)...);            
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cyclesList.push_back(cycles);
    }
    total_cycles /= controls.REP;

    cyclesList.sort();
    cycles = total_cycles;
    return  cycles;
}

template<typename comp_func,
         typename ... Args>
double measure_cycles(comp_func f, BenchControls controls, Args&&... args) {
    return measure_cycles(f, nullptr, controls, std::forward<Args>(args)...);
}

template<typename comp_func,
         typename ... Args>
double measure_cycles(comp_func f, Args&&... args) {
    BenchControls controls;
    return measure_cycles(f, nullptr, controls, std::forward<Args>(args)...);
}

template<typename comp_func>
double measure_cycles(comp_func f) {
    BenchControls controls;
    return measure_cycles(f, controls);
}

template<typename comp_func>
double measure_cycles(comp_func f, BenchControls controls) {
    double cycles = 0.;
    long num_runs = controls.NUM_RUNS;
    double multiplier = 1;
    myInt64 start, end;

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
#ifdef WARMUP
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            f();          
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (controls.CYCLES_REQUIRED) / (cycles);
        
    } while (multiplier > 2);
#endif

    std::list<double> cyclesList;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < controls.REP; j++) {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            f();            
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cyclesList.push_back(cycles);
    }
    total_cycles /= controls.REP;

    cyclesList.sort();
    cycles = total_cycles;
    return  cycles;
}

