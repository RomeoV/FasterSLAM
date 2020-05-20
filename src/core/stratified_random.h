#pragma once
#include "typedefs.h"
#include <stddef.h>

/*!
    Stratified random vector (Sampling with subpopulations).
    @param[in]  N     Interval size (inverse).
    @param[out] di    Generated array of random stratified numbers.
 */
void stratified_random(const size_t N, double* di);
double stratified_random_flops(const size_t N, double* di);
double stratified_random_memory(const size_t N, double* di);

void stratified_random_base(const size_t N, double* di);
double stratified_random_base_flops(const size_t N, double* di);
double stratified_random_base_memory(const size_t N, double* di);

/*!
    Generates a random, uniformly sampled number between [0,1].
    @return random uniform number in [0,1].
 */
double unifRand();
double unifRand_flops();
double unifRand_memory();