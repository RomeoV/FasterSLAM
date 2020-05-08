#pragma once

static unsigned long x=123456789, y=362436069, z=521288629;
unsigned long ulong_max = 18446744073709551615;

//! Pseudo RNG
//! Taken from https://stackoverflow.com/a/1640399/5616591
unsigned long xorshf96(void) { //period 2^96-1
    unsigned long t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

    t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;

    return z;
}