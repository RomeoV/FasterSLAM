#pragma once

#include <immintrin.h>
#include <cstdint>
#include <iostream>
#include <sys/time.h>



struct xorshift128plus_key_s {
    uint64_t part1;
    uint64_t part2;
};

typedef struct xorshift128plus_key_s xorshift128plus_key_t;



typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;



struct avx_xorshift128plus_key_s {
    __m256i part1;
    __m256i part2;
};

typedef struct avx_xorshift128plus_key_s avx_xorshift128plus_key_t;


typedef struct avx2_pcg_state_setseq_64 {
  __m256i state[2]; // RNG state.  All values are possible.
  __m256i inc[2];   // Controls which RNG sequence (stream) is selected. Must
                    // *always* be odd.
  __m256i
      pcg32_mult_l = _mm256_set1_epi64x((long long) (0x5851f42d4c957f2dULL & 0xffffffffu)); // set to _mm256_set1_epi64x(UINT64_C(0x5851f42d4c957f2dULL)
                    // & 0xffffffffu);
  __m256i
      pcg32_mult_h = _mm256_set1_epi64x((long long) (0x5851f42d4c957f2dULL >> 32));; // set to _mm256_set1_epi64x(UINT64_C(0x5851f42d4c957f2dULL)
                    // >> 32);

} avx2_pcg32_random_t;


// Variables

extern pcg32_random_t pcg32_seed;
extern xorshift128plus_key_t key;

extern uint64_t g_seed;
extern int rand_list[1000];

extern unsigned long long wyseed;
extern __uint128_t g_lehmer64_state;
extern uint64_t wyhash64_x; 
extern int read_rand_cntr;

extern unsigned long x,y,z;
extern unsigned long ulong_max;

extern __m256i avx_g_seed;

extern const __m256i avx_g_mult;
extern const __m256i avx_g_add;

extern const __m256i avx_uint64_max;
extern avx_xorshift128plus_key_t avx_xorshift128plus_seed;
extern avx2_pcg32_random_t avx2_pcg32_seed;

// FUNCTIONS


static inline void xorshift128plus_srand(uint64_t p1, uint64_t p2) {
  key.part1 = p1;
  key.part2 = p2;
}

static uint64_t xorshift128plus() {
    uint64_t s1 = key.part1;
    const uint64_t s0 = key.part2;
    key.part2 = s0;
    s1 ^= s1 << 23; // a
    key.part2 = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); // b, c
    return key.part2 + s0;
}

static inline void pcg32_srand(uint64_t st, uint64_t inc) {
  pcg32_seed.state=st;
  pcg32_seed.inc = inc;
}

static uint32_t pcg32()
{
    uint64_t oldstate = pcg32_seed.state;
    // Advance internal state
    pcg32_seed.state = oldstate * 6364136223846793005ULL + pcg32_seed.inc;
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// Used to seed the generator.           
static inline void fast_srand(int seed) {
    g_seed = seed;
}

// Compute a pseudorandom integer.
// Output value in range [0, 32767]
static inline int fast_rand(void) {
    g_seed = (214013*g_seed+2531011);
    return (g_seed>>32) & 0x7FFFFFFF;
}

static inline void wyrng_srand(int seed){
  wyseed = seed;
}

static inline unsigned long long wyrng(void) {
  wyseed += 0x60bee2bee120fc15ull;
  __uint128_t tmp = (__uint128_t)(wyseed)*0xa3b195354a39b70dull;
  unsigned long long m1 = (tmp >> 64) ^ tmp;
  tmp = (__uint128_t)m1 * 0x1b03738712fad5c9ull;
  return (tmp >> 64) ^ tmp;
}




static inline void lehmer64_srand(int seed){
  /* ASSERT IF SEED 0!!!*/
  g_lehmer64_state = seed;
}


static inline uint32_t lehmer64() {
  g_lehmer64_state *= 0xda942042e4dd58b5;
  return g_lehmer64_state >> 64;
}


static inline void wyhash64_srand(int seed){
  wyhash64_x = seed;
}

static inline uint64_t wyhash64() {
  wyhash64_x += 0x60bee2bee120fc15;
  __uint128_t tmp;
  tmp = (__uint128_t) wyhash64_x * 0xa3b195354a39b70d;
  uint64_t m1 = (tmp >> 64) ^ tmp;
  tmp = (__uint128_t)m1 * 0x1b03738712fad5c9;
  uint64_t m2 = (tmp >> 64) ^ tmp;
  return m2;
}



//! Pseudo RNG
//! Taken from https://stackoverflow.com/a/1640399/5616591
static inline unsigned long xorshf96(void) { //period 2^96-1
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

//////////////////// SIMD ////////////////////////////////////////

static __m256i mul64_haswell (__m256i a, __m256i b) {
    // instruction does not exist. Split into 32-bit multiplies
    __m256i bswap   = _mm256_shuffle_epi32(b,0xB1);           // swap H<->L
    __m256i prodlh  = _mm256_mullo_epi32(a,bswap);            // 32 bit L*H products

    // or use pshufb instead of psrlq to reduce port0 pressure on Haswell
    __m256i prodlh2 = _mm256_srli_epi64(prodlh, 32);          // 0  , a0Hb0L,          0, a1Hb1L
    __m256i prodlh3 = _mm256_add_epi32(prodlh2, prodlh);      // xxx, a0Lb0H+a0Hb0L, xxx, a1Lb1H+a1Hb1L
    __m256i prodlh4 = _mm256_and_si256(prodlh3, _mm256_set1_epi64x(0x00000000FFFFFFFF)); // zero high halves

    __m256i prodll  = _mm256_mul_epu32(a,b);                  // a0Lb0L,a1Lb1L, 64 bit unsigned products
    __m256i prod    = _mm256_add_epi64(prodll,prodlh4);       // a0Lb0L+(a0Lb0H+a0Hb0L)<<32, a1Lb1L+(a1Lb1H+a1Hb1L)<<32
    return  prod;
}

static inline void avx_fast_srand(uint64_t s0, uint64_t s1, uint64_t s2, uint64_t s3) {
  avx_g_seed = _mm256_set_epi64x(s3,s2,s1,s0);
}

static __m256i avx_fast_rand(void) {
  avx_g_seed = mul64_haswell(avx_g_seed, avx_g_mult);
  avx_g_seed = _mm256_add_epi64(avx_g_seed, avx_g_add);

  auto ymm0 = _mm256_srli_epi64(avx_g_seed, 32);

  return _mm256_and_si256(ymm0, avx_uint64_max);
}


/* used by xorshift128plus_jump_onkeys */
static void xorshift128plus_onkeys(uint64_t * ps0, uint64_t * ps1) {
	uint64_t s1 = *ps0;
	const uint64_t s0 = *ps1;
	*ps0 = s0;
	s1 ^= s1 << 23; // a
	*ps1 = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); // b, c
}

/* used by avx_xorshift128plus_init */
static void xorshift128plus_jump_onkeys(uint64_t in1, uint64_t in2,
		uint64_t * output1, uint64_t * output2) {
	/* todo: streamline */
	static const uint64_t JUMP[] = { 0x8a5cd789635d2dff, 0x121fd2155c472f96 };
	uint64_t s0 = 0;
	uint64_t s1 = 0;
	for (unsigned int i = 0; i < sizeof(JUMP) / sizeof(*JUMP); i++)
		for (int b = 0; b < 64; b++) {
			if (JUMP[i] & 1ULL << b) {
				s0 ^= in1;
				s1 ^= in2;
			}
			xorshift128plus_onkeys(&in1, &in2);
		}
	output1[0] = s0;
	output2[0] = s1;
}

static void avx_xorshift128plus_init(uint64_t key1, uint64_t key2) {
	uint64_t S0[4];
	uint64_t S1[4];
	S0[0] = key1;
	S1[0] = key2;
	xorshift128plus_jump_onkeys(*S0, *S1, S0 + 1, S1 + 1);
	xorshift128plus_jump_onkeys(*(S0 + 1), *(S1 + 1), S0 + 2, S1 + 2);
	xorshift128plus_jump_onkeys(*(S0 + 2), *(S1 + 2), S0 + 3, S1 + 3);
	avx_xorshift128plus_seed.part1 = _mm256_loadu_si256((const __m256i *) S0);
	avx_xorshift128plus_seed.part2 = _mm256_loadu_si256((const __m256i *) S1);
}

/*
 Return a 256-bit random "number"
 */
static __m256i avx_xorshift128plus() {
	__m256i s1 = avx_xorshift128plus_seed.part1;
	const __m256i s0 = avx_xorshift128plus_seed.part2;
	avx_xorshift128plus_seed.part1 = avx_xorshift128plus_seed.part2;
	s1 = _mm256_xor_si256(avx_xorshift128plus_seed.part2, _mm256_slli_epi64(avx_xorshift128plus_seed.part2, 23));
	avx_xorshift128plus_seed.part2 = _mm256_xor_si256(
			_mm256_xor_si256(_mm256_xor_si256(s1, s0),
					_mm256_srli_epi64(s1, 18)), _mm256_srli_epi64(s0, 5));
	return _mm256_add_epi64(avx_xorshift128plus_seed.part2, s0);
}

// credit Wenzel Jakob and Daniel Lemire 
// https://github.com/wjakob/pcg32/blob/master/pcg32_8.h and https://github.com/lemire/simdpcg/blob/master/include/simdpcg32.h
static inline __m256i avx2_pcg32() {
  const __m256i mask_l = _mm256_set1_epi64x(UINT64_C(0x00000000ffffffff));
  const __m256i shift0 = _mm256_set_epi32(7, 7, 7, 7, 6, 4, 2, 0);
  const __m256i shift1 = _mm256_set_epi32(6, 4, 2, 0, 7, 7, 7, 7);
  const __m256i const32 = _mm256_set1_epi32(32);

  __m256i s0 = avx2_pcg32_seed.state[0], s1 = avx2_pcg32_seed.state[1];

  /* Extract low and high words for partial products below */
  __m256i s0_l = _mm256_and_si256(s0, mask_l);
  __m256i s0_h = _mm256_srli_epi64(s0, 32);
  __m256i s1_l = _mm256_and_si256(s1, mask_l);
  __m256i s1_h = _mm256_srli_epi64(s1, 32);

  /* Improve high bits using xorshift step */
  __m256i s0s = _mm256_srli_epi64(s0, 18);
  __m256i s1s = _mm256_srli_epi64(s1, 18);

  __m256i s0x = _mm256_xor_si256(s0s, s0);
  __m256i s1x = _mm256_xor_si256(s1s, s1);

  __m256i s0xs = _mm256_srli_epi64(s0x, 27);
  __m256i s1xs = _mm256_srli_epi64(s1x, 27);

  __m256i xors0 = _mm256_and_si256(mask_l, s0xs);
  __m256i xors1 = _mm256_and_si256(mask_l, s1xs);

  /* Use high bits to choose a bit-level rotation */
  __m256i rot0 = _mm256_srli_epi64(s0, 59);
  __m256i rot1 = _mm256_srli_epi64(s1, 59);

  /* 64 bit multiplication using 32 bit partial products :( */
  __m256i m0_hl = _mm256_mul_epu32(s0_h, avx2_pcg32_seed.pcg32_mult_l);
  __m256i m1_hl = _mm256_mul_epu32(s1_h, avx2_pcg32_seed.pcg32_mult_l);
  __m256i m0_lh = _mm256_mul_epu32(s0_l, avx2_pcg32_seed.pcg32_mult_h);
  __m256i m1_lh = _mm256_mul_epu32(s1_l, avx2_pcg32_seed.pcg32_mult_h);

  /* Assemble lower 32 bits, will be merged into one 256 bit vector below */
  xors0 = _mm256_permutevar8x32_epi32(xors0, shift0);
  rot0 = _mm256_permutevar8x32_epi32(rot0, shift0);
  xors1 = _mm256_permutevar8x32_epi32(xors1, shift1);
  rot1 = _mm256_permutevar8x32_epi32(rot1, shift1);

  /* Continue with partial products */
  __m256i m0_ll = _mm256_mul_epu32(s0_l, avx2_pcg32_seed.pcg32_mult_l);
  __m256i m1_ll = _mm256_mul_epu32(s1_l, avx2_pcg32_seed.pcg32_mult_l);

  __m256i m0h = _mm256_add_epi64(m0_hl, m0_lh);
  __m256i m1h = _mm256_add_epi64(m1_hl, m1_lh);

  __m256i m0hs = _mm256_slli_epi64(m0h, 32);
  __m256i m1hs = _mm256_slli_epi64(m1h, 32);

  __m256i s0n = _mm256_add_epi64(m0hs, m0_ll);
  __m256i s1n = _mm256_add_epi64(m1hs, m1_ll);

  __m256i xors = _mm256_or_si256(xors0, xors1);
  __m256i rot = _mm256_or_si256(rot0, rot1);

  avx2_pcg32_seed.state[0] = _mm256_add_epi64(s0n, avx2_pcg32_seed.inc[0]);
  avx2_pcg32_seed.state[1] = _mm256_add_epi64(s1n, avx2_pcg32_seed.inc[1]);

  /* Finally, rotate and return the result */
  __m256i result =
      _mm256_or_si256(_mm256_srlv_epi32(xors, rot),
                      _mm256_sllv_epi32(xors, _mm256_sub_epi32(const32, rot)));

  return result;
}

static void avx2_pcg32_srand(const uint64_t initstate[8], const uint64_t initseq[8]) {
        const __m256i one = _mm256_set1_epi64x((long long) 1);

        avx2_pcg32_seed.state[0] = avx2_pcg32_seed.state[1] = _mm256_setzero_si256();
        avx2_pcg32_seed.inc[0] = _mm256_or_si256(
            _mm256_slli_epi64(_mm256_load_si256((__m256i *) &initseq[0]), 1),
            one);
        avx2_pcg32_seed.inc[1] = _mm256_or_si256(
            _mm256_slli_epi64(_mm256_load_si256((__m256i *) &initseq[4]), 1),
            one);
        avx2_pcg32();

        avx2_pcg32_seed.state[0] = _mm256_add_epi64(avx2_pcg32_seed.state[0], _mm256_load_si256((__m256i *) &initstate[0]));
        avx2_pcg32_seed.state[1] = _mm256_add_epi64(avx2_pcg32_seed.state[1], _mm256_load_si256((__m256i *) &initstate[4]));
        avx2_pcg32();
}

static inline int read_rand() {
  
  int temp = rand_list[read_rand_cntr];
  read_rand_cntr++;
  if (read_rand_cntr >= 1000) {
    read_rand_cntr = 0;
  }
  return temp;
}