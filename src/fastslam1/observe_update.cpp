#include "observe_update.h"

#include "add_control_noise.h"
#include "add_feature.h"
#include "add_observation_noise.h"
#include "compute_weight.h"
#include "data_associate_known.h"
#include "fastslam1_sim.h"
#include "feature_update.h"
#include "get_observations.h"
#include "linalg.h"
#include "predict.h"
#include "resample_particles.h"
#include "fastslam1_utils.h"
#include <math.h>
#include <assert.h>
#include "configfile.h"
#include "tscheb_sine.h"

#include <stddef.h>
#include <immintrin.h>


__m256d two_pi_vec = _mm256_set1_pd(2.0*M_PI);
__m256d minus_half = _mm256_set1_pd(-0.5);


__m256i z_load_mask = _mm256_set_epi64x(6,4,2,0);


//For exponential
const __m256d c0 =  _mm256_set1_pd (0.041944388f);
const __m256d c1 =  _mm256_set1_pd (0.168006673f);
const __m256d c2 =  _mm256_set1_pd (0.499999940f);
const __m256d c3 =  _mm256_set1_pd (0.999956906f);
const __m256d c4 =  _mm256_set1_pd (0.999999642f);

const __m256d l2e = _mm256_set1_pd (1.442695041); /* log2(e) */
const __m256d l2h = _mm256_set1_pd (-6.93145752e-1); /* -log(2)_hi */
const __m256d l2l = _mm256_set1_pd (-1.42860677e-6); /* -log(2)_lo */

__m256d zeros = _mm256_set1_pd(0.0);


__m256d exp_avx2_pd (__m256d x)
{
    //Source: https://stackoverflow.com/questions/48863719/fastest-implementation-of-exponential-function-using-avx

    // Philipp L.: Adapted to double-precision and added a check that the result is non-negative
    __m256d t, f, p, r;
    __m256i i, j;

    /* coefficients for core approximation to exp() in [-log(2)/2, log(2)/2] */

    // Swap to these coefficients for slighly higher accuracy
    // const __m256d c0 =  _mm256_set1_pd (0.008301110);
    // const __m256d c1 =  _mm256_set1_pd (0.041906696);
    // const __m256d c2 =  _mm256_set1_pd (0.166674897);
    // const __m256d c3 =  _mm256_set1_pd (0.499990642);
    // const __m256d c4 =  _mm256_set1_pd (0.999999762);
    // const __m256d c5 =  _mm256_set1_pd (1.000000000);

    /* exp(x) = 2^i * e^f; i = rint (log2(e) * x), f = x - log(2) * i */
    t = _mm256_mul_pd (x, l2e);      /* t = log2(e) * x */
    r = _mm256_round_pd (t, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC); /* r = rint (t) */

    f = _mm256_fmadd_pd (r, l2h, x); /* x - log(2)_hi * r */
    f = _mm256_fmadd_pd (r, l2l, f); /* f = x - log(2)_hi * r - log(2)_lo * r */

    //i = _mm256_cvtpd_epi32(t);       /* i = (int)rint(t) */
    i = _mm256_cvtepi32_epi64(_mm256_cvtpd_epi32(t));

    /* p ~= exp (f), -log(2)/2 <= f <= log(2)/2 */
    p = c0;                          /* c0 */
    p = _mm256_fmadd_pd (p, f, c1);  /* c0*f+c1 */
    p = _mm256_fmadd_pd (p, f, c2);  /* (c0*f+c1)*f+c2 */
    p = _mm256_fmadd_pd (p, f, c3);  /* ((c0*f+c1)*f+c2)*f+c3 */
    p = _mm256_fmadd_pd (p, f, c4);  /* (((c0*f+c1)*f+c2)*f+c3)*f+c4 ~= exp(f) */
    // p = _mm256_fmadd_pd (p, f, c5); // Swap to these coefficients for slighly higher accuracy

    /* exp(x) = 2^i * p */
    j = _mm256_slli_epi64 (i, 52); /* i << 23 (on epi32, <<52 on epi64 (Philipp L., tested ^^)*/
    r = _mm256_castsi256_pd (_mm256_add_epi64 (j, _mm256_castpd_si256 (p))); /* r = p * 2^i */

    //Philipp L. : Function breaks at super negative exponents -> yields negative numbers. This check prevents this.
    r = _mm256_blendv_pd( r, zeros,  r);

    return r;
}

extern __inline __m256d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm256_load2_m128d (double const *__PH, double const *__PL)
{
  return _mm256_insertf128_pd (_mm256_castpd128_pd256 (_mm_loadu_pd (__PL)),
			       _mm_loadu_pd (__PH), 1);
}

void observe_update(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    observe_update_fast(lm, N_features, xtrue, R, ftag, da_table, ftag_visible, z, Nf_visible, zf, idf,
                        zn, particles, weights);
}

void observe_update_base(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    // Compute true data, then add noise
    // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
    // memcpy(ftag_visible, ftag, N_features*sizeof(int));
    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }

    //z is the range and bearing of the observed landmark
    
    get_observations_base(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
    //print(*z,*Nf_visible,2,std::cout);
    if ( *Nf_visible == 0 ) {
        return;
    }
    
    add_observation_noise_base(z, *Nf_visible, R, SWITCH_SENSOR_NOISE);

    //Compute (known) data associations
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known_base(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions

    // perform update
    for (size_t i = 0; i < NPARTICLES; i++) {
        if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
            // std::cout<<"HiBase"<<std::endl;
            // TODO: precompute the jacobian of compute_weight_base and feature_update_base
            // this will half the processing need of compute_jacobians
            Vector2d zp[count_zf] __attribute__ ((aligned(32)));
            Matrix23d Hv[count_zf] __attribute__ ((aligned(32)));
            Matrix2d Hf[count_zf] __attribute__ ((aligned(32)));
            Matrix2d Sf[count_zf] __attribute__ ((aligned(32)));
            double w = compute_weight_base(&particles[i], zf, count_zf, idf, R, zp, Hv, Hf, Sf);
            w *= weights[i];
            weights[i] = w;
            feature_update_base(&particles[i], zf, idf, count_zf, R, zp, Hv, Hf, Sf);
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature_base(&particles[i], zn, count_zn, R);
        }
    }

    resample_particles_base(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}

void observe_update_active(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    // Compute true data, then add noise
    // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
    // memcpy(ftag_visible, ftag, N_features*sizeof(int));
    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }

    //z is the range and bearing of the observed landmark
    
    get_observations(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
    //print(*z,*Nf_visible,2,std::cout);
    if ( *Nf_visible == 0 ) {
        return;
    }
    
    add_observation_noise(z, *Nf_visible, R, SWITCH_SENSOR_NOISE);

    //Compute (known) data associations
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions

    // perform update
    for (size_t i = 0; i < NPARTICLES; i++) {
        if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
            Vector2d zp[count_zf] __attribute__ ((aligned(32)));
            Matrix23d Hv[count_zf] __attribute__ ((aligned(32)));
            Matrix2d Hf[count_zf] __attribute__ ((aligned(32)));
            Matrix2d Sf[count_zf] __attribute__ ((aligned(32)));
            compute_jacobians_fast(&particles[i], idf, count_zf, R, zp, Hv, Hf, Sf);
            double w = compute_weight_active(&particles[i], zf, count_zf, idf, R, zp, Hv, Hf, Sf);
            w *= weights[i];
            weights[i] = w;
            feature_update_active(&particles[i], zf, idf, count_zf, R, zp, Hv, Hf, Sf);
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature(&particles[i], zn, count_zn, R);
        }
    }

    resample_particles(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}

void observe_update_fast(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    // Compute true data, then add noise
    // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
    // memcpy(ftag_visible, ftag, N_features*sizeof(int));
    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }
    //z is the range and bearing of the observed landmark
    get_observations(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
    //print(*z,*Nf_visible,2,std::cout);
    if ( *Nf_visible == 0 ) {
        return;
    }
    add_observation_noise(z, *Nf_visible, R, SWITCH_SENSOR_NOISE);

    //Compute (known) data associations
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions
    
    //Variables
    Matrix23d Hv[4* count_zf] __attribute__ ((aligned(32))); //Unused in whole code
    Matrix2d Hf[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf1[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf2[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf3[count_zf] __attribute__ ((aligned(32)));
    Vector2d feat_diff[4*count_zf] __attribute__ ((aligned(32)));
    Vector2d zp[4*count_zf] __attribute__ ((aligned(32)));
    Matrix2d Sf[4*count_zf] __attribute__ ((aligned(32)));

    Vector2d* zp1 = zp+ 1*count_zf;
    Vector2d* zp2 = zp+ 2*count_zf;
    Vector2d* zp3 = zp+ 3*count_zf;

    Vector2d* feat_diff1 = feat_diff+ 1*count_zf;
    Vector2d* feat_diff2 = feat_diff+ 2*count_zf;
    Vector2d* feat_diff3 = feat_diff+ 3*count_zf;

    Matrix2d* Sf1 = Sf + 1*count_zf;
    Matrix2d* Sf2 = Sf + 2 *count_zf;
    Matrix2d* Sf3 = Sf + 3*count_zf;

    double den[4]  __attribute__ ((aligned(32))), num[4]  __attribute__ ((aligned(32)));
    double vT_S_inv_v[4]  __attribute__ ((aligned(32)));

    __m256d den_v, num_v, vec_vSv, w_vec, zp_dist, zp_angle, z_dist, z_angle, fdiff_dist, fdiff_angle, fdiff1, fdiff2;
    __m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7, ymm8, ymm9, sinv0, sinv1, sinv2, sinv3, s_zeros, s_ones, s_twos, s_threes;
    __m256d determinants, inv_det, neg_inv_det, weights_v;
    __m256i zp_load_mask = _mm256_set_epi64x(2*3*count_zf,2*2*count_zf,2*1*count_zf,0); ///!!!!! Order here
    __m256i sf_mask = _mm256_set_epi64x(4*3*count_zf,4*2*count_zf,4*1*count_zf,0);

    __m256d R_vec = _mm256_load_pd(R);
    
    for (size_t i = 0; i < NPARTICLES; i+=4) {
        if ( count_zf != 0 ) { 
            //COMPUTE JACOBIANS
            compute_jacobians_fast(particles + i, idf, count_zf, R, zp, Hv, Hf, Sf);
            compute_jacobians_fast(particles + i +1, idf, count_zf, R, zp1, Hv + 1* count_zf, Hf1, Sf1);
            compute_jacobians_fast(particles + i +2, idf, count_zf, R, zp2, Hv + 2* count_zf, Hf2, Sf2);
            compute_jacobians_fast(particles + i +3, idf, count_zf, R, zp3, Hv + 3* count_zf, Hf3, Sf3);
            //END COMPUTE_JACOBIANS

            weights_v = _mm256_load_pd(weights+i);
            for (size_t j = 0; j < count_zf; j++) {
                // Particle 0
                z_dist = _mm256_set1_pd(zf[j][0]);
                z_angle = _mm256_set1_pd(zf[j][1]);

                zp_dist = _mm256_i64gather_pd(*(zp+j), zp_load_mask,8);
                zp_angle = _mm256_i64gather_pd(*(zp+j)+1, zp_load_mask,8);

                fdiff_dist = _mm256_sub_pd(z_dist, zp_dist);
                fdiff_angle = _mm256_sub_pd(z_angle, zp_angle);

                fdiff_angle = simple_pi_to_pi_avx(fdiff_angle);
                
                fdiff1 = _mm256_shuffle_pd(fdiff_dist, fdiff_angle, 0b0000);
                fdiff2 = _mm256_shuffle_pd(fdiff_dist, fdiff_angle, 0b1111);

                _mm256_store2_m128d(feat_diff2[j], feat_diff[j], fdiff1);
                _mm256_store2_m128d(feat_diff3[j], feat_diff1[j], fdiff2);

                //INV and determinant simultaneously
                s_zeros = _mm256_i64gather_pd(Sf[j], sf_mask,8); //Sfp0_0, Sf1_0 ..: Sf_0, Sf1_0, Sf2_0, Sf3_0
                s_ones = _mm256_i64gather_pd(Sf[j]+1, sf_mask,8); //Sf1_0, Sf1_1
                s_twos = _mm256_i64gather_pd(Sf[j]+2, sf_mask,8);
                s_threes = _mm256_i64gather_pd(Sf[j]+3, sf_mask,8);

                ymm0 = _mm256_mul_pd(s_zeros, s_threes);
                ymm1 = _mm256_mul_pd(s_twos, s_ones);
                
                determinants = _mm256_sub_pd(ymm0, ymm1); //dets0, dets1, dets2, dets3

                inv_det = _mm256_div_pd(_mm256_set1_pd(1.0), determinants);
                neg_inv_det = _mm256_mul_pd(_mm256_set1_pd(-1.0), inv_det);

                s_threes  = _mm256_mul_pd(inv_det, s_threes); //s * A[3]
                s_ones = _mm256_mul_pd(neg_inv_det, s_ones);
                s_twos = _mm256_mul_pd(neg_inv_det, s_twos);
                s_zeros = _mm256_mul_pd(inv_det, s_zeros);

                ymm6 = _mm256_unpacklo_pd(s_threes , s_ones);
                ymm7 = _mm256_unpackhi_pd(s_threes , s_ones);
                ymm8 = _mm256_unpacklo_pd(s_twos, s_zeros);
                ymm9 = _mm256_unpackhi_pd(s_twos, s_zeros);

                sinv0 = _mm256_permute2f128_pd(ymm6, ymm8, 0b00100000);
                sinv1 = _mm256_permute2f128_pd(ymm7, ymm9, 0b00100000);
                sinv2 = _mm256_permute2f128_pd(ymm6, ymm8, 0b00110001);
                sinv3 = _mm256_permute2f128_pd(ymm7, ymm9, 0b00110001);

                //vTMv
                ymm0 = _mm256_permute2f128_pd(sinv0, sinv2, 0b00100000);
                ymm2 = _mm256_permute2f128_pd(sinv0, sinv2, 0b00110001);

                ymm1 = _mm256_permute2f128_pd(sinv1, sinv3, 0b00100000);
                ymm3 = _mm256_permute2f128_pd(sinv1, sinv3, 0b00110001);

                ymm0 = _mm256_mul_pd(fdiff1, ymm0);
                ymm2 = _mm256_mul_pd(fdiff1, ymm2);
                
                ymm1 = _mm256_mul_pd(fdiff2, ymm1);
                ymm3 = _mm256_mul_pd(fdiff2, ymm3);

                ymm0 = _mm256_hadd_pd(ymm0, ymm2);
                ymm1 = _mm256_hadd_pd(ymm1, ymm3);

                ymm2 = _mm256_mul_pd(fdiff1, ymm0);
                ymm3 = _mm256_mul_pd(fdiff2, ymm1);
                vec_vSv =  _mm256_hadd_pd(ymm2, ymm3);
                //vTMv
                
                //Update Weights
                ymm0 = _mm256_sqrt_pd(determinants);
                den_v = _mm256_mul_pd(ymm0, two_pi_vec);
                num_v = exp_avx2_pd(_mm256_mul_pd(minus_half, vec_vSv));

                weights_v = _mm256_mul_pd(weights_v, _mm256_div_pd(num_v, den_v));

#ifdef KF_YGLEE
                KF_cholesky_update(particles[i].xf + 2 * idf[j], particles[i].Pf + 4 * idf[j], 
                                feat_diff[j], R, 
                                Hf[j]);
                KF_cholesky_update(particles[i+1].xf + 2 * idf[j], particles[i+1].Pf + 4 * idf[j], 
                                feat_diff1[j], R, 
                                Hf1[j]);
                KF_cholesky_update(particles[i+2].xf + 2 * idf[j], particles[i+2].Pf + 4 * idf[j], 
                                feat_diff2[j], R, 
                                Hf2[j]);
                KF_cholesky_update(particles[i+3].xf + 2 * idf[j], particles[i+3].Pf + 4 * idf[j], 
                                feat_diff3[j], R, 
                                Hf3[j]);
#else
                __m256d xfp0p2 = _mm256_load2_m128d(particles[i+2].xf + 2 * idf[j], particles[i].xf + 2 * idf[j]);
                __m256d xfp1p3 = _mm256_load2_m128d(particles[i+3].xf + 2 * idf[j], particles[i+1].xf + 2 * idf[j]);

                __m256d Pf0 = _mm256_load_pd(particles[i].Pf + 4 * idf[j]);
                __m256d Pf1 = _mm256_load_pd(particles[i+1].Pf + 4 * idf[j]);
                __m256d Pf2 = _mm256_load_pd(particles[i+2].Pf + 4 * idf[j]);
                __m256d Pf3 = _mm256_load_pd(particles[i+3].Pf + 4 * idf[j]);

                __m256d Hf0v = _mm256_load_pd(Hf[j]);
                __m256d Hf1v = _mm256_load_pd(Hf1[j]);
                __m256d Hf2v = _mm256_load_pd(Hf2[j]);
                __m256d Hf3v = _mm256_load_pd(Hf3[j]);

                KF_cholesky_update_unrolled4_avx(&xfp0p2, &xfp1p3, &Pf0, &Pf1, &Pf2, &Pf3, 
                                fdiff1, fdiff2, R_vec, Hf0v, Hf1v, Hf2v, Hf3v);
                
                _mm256_store2_m128d(particles[i+2].xf + 2 * idf[j], particles[i].xf + 2 * idf[j], xfp0p2);
                _mm256_store2_m128d(particles[i+3].xf + 2 * idf[j], particles[i+1].xf + 2 * idf[j], xfp1p3);

                _mm256_store_pd(particles[i].Pf + 4 * idf[j], Pf0);
                _mm256_store_pd(particles[i+1].Pf + 4 * idf[j], Pf1);
                _mm256_store_pd(particles[i+2].Pf + 4 * idf[j], Pf2);
                _mm256_store_pd(particles[i+3].Pf + 4 * idf[j], Pf3);     
#endif
            }
            _mm256_store_pd(weights+i, weights_v); 
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature(particles+i, zn, count_zn, R);
            add_feature(particles+i+1, zn, count_zn, R);
            add_feature(particles+i+2, zn, count_zn, R);
            add_feature(particles+i+3, zn, count_zn, R);
        }
    }
    resample_particles(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}

void observe_update_fast_KF_comp_not_unrolled(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    // Compute true data, then add noise
    // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
    // memcpy(ftag_visible, ftag, N_features*sizeof(int));
    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }
    //z is the range and bearing of the observed landmark
    get_observations(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
    //print(*z,*Nf_visible,2,std::cout);
    if ( *Nf_visible == 0 ) {
        return;
    }
    add_observation_noise(z, *Nf_visible, R, SWITCH_SENSOR_NOISE);

    //Compute (known) data associations
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions
    
    //Variables
    Matrix23d Hv[4* count_zf] __attribute__ ((aligned(32))); //Unused in whole code
    Matrix2d Hf[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf1[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf2[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf3[count_zf] __attribute__ ((aligned(32)));
    Vector2d feat_diff[4*count_zf] __attribute__ ((aligned(32)));
    Vector2d zp[4*count_zf] __attribute__ ((aligned(32)));
    Matrix2d Sf[4*count_zf] __attribute__ ((aligned(32)));

    Vector2d* zp1 = zp+ 1*count_zf;
    Vector2d* zp2 = zp+ 2*count_zf;
    Vector2d* zp3 = zp+ 3*count_zf;

    Vector2d* feat_diff1 = feat_diff+ 1*count_zf;
    Vector2d* feat_diff2 = feat_diff+ 2*count_zf;
    Vector2d* feat_diff3 = feat_diff+ 3*count_zf;

    Matrix2d* Sf1 = Sf + 1*count_zf;
    Matrix2d* Sf2 = Sf + 2 *count_zf;
    Matrix2d* Sf3 = Sf + 3*count_zf;

    double den[4]  __attribute__ ((aligned(32))), num[4]  __attribute__ ((aligned(32)));
    double vT_S_inv_v[4]  __attribute__ ((aligned(32)));

    __m256d den_v, num_v, vec_vSv, w_vec, zp_dist, zp_angle, z_dist, z_angle, fdiff_dist, fdiff_angle, fdiff1, fdiff2;
    __m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7, ymm8, ymm9, sinv0, sinv1, sinv2, sinv3, s_zeros, s_ones, s_twos, s_threes;
    __m256d determinants, inv_det, neg_inv_det, weights_v;
    __m256i zp_load_mask = _mm256_set_epi64x(2*3*count_zf,2*2*count_zf,2*1*count_zf,0); ///!!!!! Order here
    __m256i sf_mask = _mm256_set_epi64x(4*3*count_zf,4*2*count_zf,4*1*count_zf,0);
    
    for (size_t i = 0; i < NPARTICLES; i+=4) {
        if ( count_zf != 0 ) { 
            //COMPUTE JACOBIANS
            compute_jacobians_fast(particles + i, idf, count_zf, R, zp, Hv, Hf, Sf);
            compute_jacobians_fast(particles + i +1, idf, count_zf, R, zp1, Hv + 1* count_zf, Hf1, Sf1);
            compute_jacobians_fast(particles + i +2, idf, count_zf, R, zp2, Hv + 2* count_zf, Hf2, Sf2);
            compute_jacobians_fast(particles + i +3, idf, count_zf, R, zp3, Hv + 3* count_zf, Hf3, Sf3);
            //END COMPUTE_JACOBIANS

            weights_v = _mm256_load_pd(weights+i);
            for (size_t j = 0; j < count_zf; j++) {
                // Particle 0
                z_dist = _mm256_set1_pd(zf[j][0]);
                z_angle = _mm256_set1_pd(zf[j][1]);

                zp_dist = _mm256_i64gather_pd(*(zp+j), zp_load_mask,8);
                zp_angle = _mm256_i64gather_pd(*(zp+j)+1, zp_load_mask,8);

                fdiff_dist = _mm256_sub_pd(z_dist, zp_dist);
                fdiff_angle = _mm256_sub_pd(z_angle, zp_angle);

                fdiff_angle = simple_pi_to_pi_avx(fdiff_angle);
                
                fdiff1 = _mm256_shuffle_pd(fdiff_dist, fdiff_angle, 0b0000);
                fdiff2 = _mm256_shuffle_pd(fdiff_dist, fdiff_angle, 0b1111);

                _mm256_store2_m128d(feat_diff2[j], feat_diff[j], fdiff1);
                _mm256_store2_m128d(feat_diff3[j], feat_diff1[j], fdiff2);

                //INV and determinant simultaneously
                s_zeros = _mm256_i64gather_pd(Sf[j], sf_mask,8); //Sfp0_0, Sf1_0 ..: Sf_0, Sf1_0, Sf2_0, Sf3_0
                s_ones = _mm256_i64gather_pd(Sf[j]+1, sf_mask,8); //Sf1_0, Sf1_1
                s_twos = _mm256_i64gather_pd(Sf[j]+2, sf_mask,8);
                s_threes = _mm256_i64gather_pd(Sf[j]+3, sf_mask,8);

                ymm0 = _mm256_mul_pd(s_zeros, s_threes);
                ymm1 = _mm256_mul_pd(s_twos, s_ones);
                
                determinants = _mm256_sub_pd(ymm0, ymm1); //dets0, dets1, dets2, dets3

                inv_det = _mm256_div_pd(_mm256_set1_pd(1.0), determinants);
                neg_inv_det = _mm256_mul_pd(_mm256_set1_pd(-1.0), inv_det);

                s_threes  = _mm256_mul_pd(inv_det, s_threes); //s * A[3]
                s_ones = _mm256_mul_pd(neg_inv_det, s_ones);
                s_twos = _mm256_mul_pd(neg_inv_det, s_twos);
                s_zeros = _mm256_mul_pd(inv_det, s_zeros);

                ymm6 = _mm256_unpacklo_pd(s_threes , s_ones);
                ymm7 = _mm256_unpackhi_pd(s_threes , s_ones);
                ymm8 = _mm256_unpacklo_pd(s_twos, s_zeros);
                ymm9 = _mm256_unpackhi_pd(s_twos, s_zeros);

                sinv0 = _mm256_permute2f128_pd(ymm6, ymm8, 0b00100000);
                sinv1 = _mm256_permute2f128_pd(ymm7, ymm9, 0b00100000);
                sinv2 = _mm256_permute2f128_pd(ymm6, ymm8, 0b00110001);
                sinv3 = _mm256_permute2f128_pd(ymm7, ymm9, 0b00110001);

                //vTMv
                ymm0 = _mm256_permute2f128_pd(sinv0, sinv2, 0b00100000);
                ymm2 = _mm256_permute2f128_pd(sinv0, sinv2, 0b00110001);

                ymm1 = _mm256_permute2f128_pd(sinv1, sinv3, 0b00100000);
                ymm3 = _mm256_permute2f128_pd(sinv1, sinv3, 0b00110001);

                ymm0 = _mm256_mul_pd(fdiff1, ymm0);
                ymm2 = _mm256_mul_pd(fdiff1, ymm2);
                
                ymm1 = _mm256_mul_pd(fdiff2, ymm1);
                ymm3 = _mm256_mul_pd(fdiff2, ymm3);

                ymm0 = _mm256_hadd_pd(ymm0, ymm2);
                ymm1 = _mm256_hadd_pd(ymm1, ymm3);

                ymm2 = _mm256_mul_pd(fdiff1, ymm0);
                ymm3 = _mm256_mul_pd(fdiff2, ymm1);
                vec_vSv =  _mm256_hadd_pd(ymm2, ymm3);
                //vTMv
                
                //Update Weights
                ymm0 = _mm256_sqrt_pd(determinants);
                den_v = _mm256_mul_pd(ymm0, two_pi_vec);
                num_v = exp_avx2_pd(_mm256_mul_pd(minus_half, vec_vSv));

                weights_v = _mm256_mul_pd(weights_v, _mm256_div_pd(num_v, den_v));

                KF_cholesky_update(particles[i].xf + 2 * idf[j], particles[i].Pf + 4 * idf[j], 
                                feat_diff[j], R, 
                                Hf[j]);
                KF_cholesky_update(particles[i+1].xf + 2 * idf[j], particles[i+1].Pf + 4 * idf[j], 
                                feat_diff1[j], R, 
                                Hf1[j]);
                KF_cholesky_update(particles[i+2].xf + 2 * idf[j], particles[i+2].Pf + 4 * idf[j], 
                                feat_diff2[j], R, 
                                Hf2[j]);
                KF_cholesky_update(particles[i+3].xf + 2 * idf[j], particles[i+3].Pf + 4 * idf[j], 
                                feat_diff3[j], R, 
                                Hf3[j]);
            }
            _mm256_store_pd(weights+i, weights_v); 
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature(particles+i, zn, count_zn, R);
            add_feature(particles+i+1, zn, count_zn, R);
            add_feature(particles+i+2, zn, count_zn, R);
            add_feature(particles+i+3, zn, count_zn, R);
        }
    }
    resample_particles(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}


void observe_update_fast_romeo_vTMv(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    // Compute true data, then add noise
    // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
    // memcpy(ftag_visible, ftag, N_features*sizeof(int));
    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }

    //z is the range and bearing of the observed landmark
    
    get_observations(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
    //print(*z,*Nf_visible,2,std::cout);
    if ( *Nf_visible == 0 ) {
        return;
    }
    
    add_observation_noise(z, *Nf_visible, R, SWITCH_SENSOR_NOISE);

    //Compute (known) data associations
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions
    //Vector2d zp[count_zf] __attribute__ ((aligned(32)));
    Matrix23d Hv[4* count_zf] __attribute__ ((aligned(32))); //Unused in whole code
    Matrix2d Hf[count_zf] __attribute__ ((aligned(32)));
    
    

    //Vector2d zp1[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf1[count_zf] __attribute__ ((aligned(32)));
    // Matrix2d Sf1[count_zf] __attribute__ ((aligned(32)));
    //Vector2d feat_diff1[count_zf] __attribute__ ((aligned(32)));

    // Vector2d zp2[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf2[count_zf] __attribute__ ((aligned(32)));
    // Matrix2d Sf2[count_zf] __attribute__ ((aligned(32)));
    //Vector2d feat_diff2[count_zf] __attribute__ ((aligned(32)));

    // Vector2d zp3[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf3[count_zf] __attribute__ ((aligned(32)));
    // Matrix2d Sf3[count_zf] __attribute__ ((aligned(32)));
    // Vector2d feat_diff3[count_zf] __attribute__ ((aligned(32)));

    Vector2d feat_diff[4*count_zf] __attribute__ ((aligned(32)));
    Vector2d zp[4*count_zf] __attribute__ ((aligned(32)));
    Matrix2d Sf[4*count_zf] __attribute__ ((aligned(32)));

    Vector2d* zp1 = zp+ 1*count_zf;
    Vector2d* zp2 = zp+ 2*count_zf;
    Vector2d* zp3 = zp+ 3*count_zf;

    Vector2d* feat_diff1 = feat_diff+ 1*count_zf;
    Vector2d* feat_diff2 = feat_diff+ 2*count_zf;
    Vector2d* feat_diff3 = feat_diff+ 3*count_zf;

    Matrix2d* Sf1 = Sf + 1*count_zf;
    Matrix2d* Sf2 = Sf + 2 *count_zf;
    Matrix2d* Sf3 = Sf + 3*count_zf;

    Vector2d S_inv_v_al[4] __attribute__ ((aligned(32))); //Order: 0,2,1,3 !!!
    

    //double dx, dy, d2, d, dinv, d2inv, dx_d2inv, dy_d2inv, dx_dinv, dy_dinv;
    //double den, num;
    Matrix2d S_inv __attribute__ ((aligned(32)));
    Vector2d S_inv_v  __attribute__ ((aligned(32)));

    Matrix2d S_inv1 __attribute__ ((aligned(32)));
    Vector2d S_inv_v1  __attribute__ ((aligned(32)));

    Matrix2d S_inv2 __attribute__ ((aligned(32)));
    Vector2d S_inv_v2  __attribute__ ((aligned(32)));

    Matrix2d S_inv3 __attribute__ ((aligned(32)));
    Vector2d S_inv_v3  __attribute__ ((aligned(32)));


    double den[4]  __attribute__ ((aligned(32))), num[4]  __attribute__ ((aligned(32)));

    double vT_S_inv_v[4]  __attribute__ ((aligned(32)));

    __m256d den_v, num_v, vec_vSv, w_vec, zp_dist, zp_angle, z_dist, z_angle, fdiff_dist, fdiff_angle, fdiff1, fdiff2;

    __m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7, ymm8, ymm9, sinv0, sinv1, sinv2, sinv3, s_zeros, s_ones, s_twos, s_threes;

    __m256d determinants, inv_det, neg_inv_det, weights_v;

    __m256i zp_load_mask = _mm256_set_epi64x(2*3*count_zf,2*2*count_zf,2*1*count_zf,0); ///!!!!! Order here
    __m256i sf_mask = _mm256_set_epi64x(4*3*count_zf,4*2*count_zf,4*1*count_zf,0);
    
    for (size_t i = 0; i < NPARTICLES; i+=4) {
        if ( count_zf != 0 ) { 
            //COMPUTE JACOBIANS
            compute_jacobians_fast(particles + i, idf, count_zf, R, zp, Hv, Hf, Sf);
            compute_jacobians_fast(particles + i +1, idf, count_zf, R, zp1, Hv + 1* count_zf, Hf1, Sf1);
            compute_jacobians_fast(particles + i +2, idf, count_zf, R, zp2, Hv + 2* count_zf, Hf2, Sf2);
            compute_jacobians_fast(particles + i +3, idf, count_zf, R, zp3, Hv + 3* count_zf, Hf3, Sf3);
            //END COMPUTE_JACOBIANS

            weights_v = _mm256_load_pd(weights+i);

            for (size_t j = 0; j < count_zf; j++) {
                // Particle 0
                z_dist = _mm256_set1_pd(zf[j][0]);
                z_angle = _mm256_set1_pd(zf[j][1]);

                zp_dist = _mm256_i64gather_pd(*(zp+j), zp_load_mask,8);
                zp_angle = _mm256_i64gather_pd(*(zp+j)+1, zp_load_mask,8);

                fdiff_dist = _mm256_sub_pd(z_dist, zp_dist);
                fdiff_angle = _mm256_sub_pd(z_angle, zp_angle);

                fdiff_angle = simple_pi_to_pi_avx(fdiff_angle);
                
                fdiff1 = _mm256_shuffle_pd(fdiff_dist, fdiff_angle, 0b0000);
                fdiff2 = _mm256_shuffle_pd(fdiff_dist, fdiff_angle, 0b1111);

                _mm256_store2_m128d(feat_diff2[j], feat_diff[j], fdiff1);
                _mm256_store2_m128d(feat_diff3[j], feat_diff1[j], fdiff2);

                //INV and determinant simultaneously
                s_zeros = _mm256_i64gather_pd(Sf[j], sf_mask,8); //Sfp0_0, Sf1_0 ..: Sf_0, Sf1_0, Sf2_0, Sf3_0
                s_ones = _mm256_i64gather_pd(Sf[j]+1, sf_mask,8); //Sf1_0, Sf1_1
                s_twos = _mm256_i64gather_pd(Sf[j]+2, sf_mask,8);
                s_threes = _mm256_i64gather_pd(Sf[j]+3, sf_mask,8);



                ymm0 = _mm256_mul_pd(s_zeros, s_threes);
                ymm1 = _mm256_mul_pd(s_twos, s_ones);
                
                determinants = _mm256_sub_pd(ymm0, ymm1); //dets0, dets1, dets2, dets3

                inv_det = _mm256_div_pd(_mm256_set1_pd(1.0), determinants);
                neg_inv_det = _mm256_mul_pd(_mm256_set1_pd(-1.0), inv_det);

                s_threes  = _mm256_mul_pd(inv_det, s_threes); //s * A[3]
                s_ones = _mm256_mul_pd(neg_inv_det, s_ones);
                s_twos = _mm256_mul_pd(neg_inv_det, s_twos);
                s_zeros = _mm256_mul_pd(inv_det, s_zeros);

                ymm6 = _mm256_unpacklo_pd(s_threes , s_ones);
                ymm7 = _mm256_unpackhi_pd(s_threes , s_ones);
                ymm8 = _mm256_unpacklo_pd(s_twos, s_zeros);
                ymm9 = _mm256_unpackhi_pd(s_twos, s_zeros);

                // _mm256_store_pd(S_inv, _mm256_permute2f128_pd(ymm6, ymm8, 0b00100000));
                // _mm256_store_pd(S_inv2,_mm256_permute2f128_pd(ymm6, ymm8, 0b00110001));
                // _mm256_store_pd(S_inv1, _mm256_permute2f128_pd(ymm7, ymm9, 0b00100000));
                // _mm256_store_pd(S_inv3, _mm256_permute2f128_pd(ymm7, ymm9, 0b00110001));

                sinv0 = _mm256_permute2f128_pd(ymm6, ymm8, 0b00100000);
                sinv1 = _mm256_permute2f128_pd(ymm7, ymm9, 0b00100000);
                sinv2 = _mm256_permute2f128_pd(ymm6, ymm8, 0b00110001);
                sinv3 = _mm256_permute2f128_pd(ymm7, ymm9, 0b00110001);


                //vTMv

                vec_vSv=mm_vT_M_v_avx2(sinv0,  sinv2,
                       sinv1, sinv3,
                       fdiff1, fdiff2);
                //vTMv

                
                //Weights
                ymm0 = _mm256_sqrt_pd(determinants);
                den_v = _mm256_mul_pd(ymm0, two_pi_vec);
                num_v = exp_avx2_pd(_mm256_mul_pd(minus_half, vec_vSv));

                weights_v = _mm256_mul_pd(weights_v, _mm256_div_pd(num_v, den_v));

                //KF Cholesky

                // KF_cholesky_update(__m256d* xfp0p2, __m256d* xfp1p3, __m256d* Pf0, __m256d* Pf1, __m256d* Pf2, __m256d* Pf3,
                //                    fdiffp0p2, fdiffp1p3, __m256d R_vec, __m256d Hfp0, __m256d Hfp1,__m256d Hfp2,__m256d Hfp3) 

                KF_cholesky_update(particles[i].xf + 2 * idf[j], particles[i].Pf + 4 * idf[j], 
                                feat_diff[j], R, 
                                Hf[j]);
                KF_cholesky_update(particles[i+1].xf + 2 * idf[j], particles[i+1].Pf + 4 * idf[j], 
                                feat_diff1[j], R, 
                                Hf1[j]);
                KF_cholesky_update(particles[i+2].xf + 2 * idf[j], particles[i+2].Pf + 4 * idf[j], 
                                feat_diff2[j], R, 
                                Hf2[j]);
                KF_cholesky_update(particles[i+3].xf + 2 * idf[j], particles[i+3].Pf + 4 * idf[j], 
                                feat_diff3[j], R, 
                                Hf3[j]);

                
            }
            _mm256_store_pd(weights+i, weights_v);
                
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature(particles+i, zn, count_zn, R);
            add_feature(particles+i+1, zn, count_zn, R);
            add_feature(particles+i+2, zn, count_zn, R);
            add_feature(particles+i+3, zn, count_zn, R);
        }
    }

    resample_particles(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}
void observe_update_fast_v1(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    // Compute true data, then add noise
    // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
    // memcpy(ftag_visible, ftag, N_features*sizeof(int));
    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }

    //z is the range and bearing of the observed landmark
    
    get_observations(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
    //print(*z,*Nf_visible,2,std::cout);
    if ( *Nf_visible == 0 ) {
        return;
    }
    
    add_observation_noise(z, *Nf_visible, R, SWITCH_SENSOR_NOISE);

    //Compute (known) data associations
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions
    
    
    //Vector2d zp[count_zf] __attribute__ ((aligned(32)));
    Matrix23d Hv[4* count_zf] __attribute__ ((aligned(32))); //Unused in whole code
    Matrix2d Hf[count_zf] __attribute__ ((aligned(32)));
    
    

    //Vector2d zp1[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf1[count_zf] __attribute__ ((aligned(32)));
    // Matrix2d Sf1[count_zf] __attribute__ ((aligned(32)));
    //Vector2d feat_diff1[count_zf] __attribute__ ((aligned(32)));

    // Vector2d zp2[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf2[count_zf] __attribute__ ((aligned(32)));
    // Matrix2d Sf2[count_zf] __attribute__ ((aligned(32)));
    //Vector2d feat_diff2[count_zf] __attribute__ ((aligned(32)));

    // Vector2d zp3[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf3[count_zf] __attribute__ ((aligned(32)));
    // Matrix2d Sf3[count_zf] __attribute__ ((aligned(32)));
    // Vector2d feat_diff3[count_zf] __attribute__ ((aligned(32)));

    Vector2d feat_diff[4*count_zf] __attribute__ ((aligned(32)));
    Vector2d zp[4*count_zf] __attribute__ ((aligned(32)));
    Matrix2d Sf[4*count_zf] __attribute__ ((aligned(32)));

    Vector2d* zp1 = zp+ 1*count_zf;
    Vector2d* zp2 = zp+ 2*count_zf;
    Vector2d* zp3 = zp+ 3*count_zf;

    Vector2d* feat_diff1 = feat_diff+ 1*count_zf;
    Vector2d* feat_diff2 = feat_diff+ 2*count_zf;
    Vector2d* feat_diff3 = feat_diff+ 3*count_zf;

    Matrix2d* Sf1 = Sf + 1*count_zf;
    Matrix2d* Sf2 = Sf + 2 *count_zf;
    Matrix2d* Sf3 = Sf + 3*count_zf;

    //double dx, dy, d2, d, dinv, d2inv, dx_d2inv, dy_d2inv, dx_dinv, dy_dinv;
    //double den, num;
    Matrix2d S_inv __attribute__ ((aligned(32)));
    Vector2d S_inv_v  __attribute__ ((aligned(32)));

    Matrix2d S_inv1 __attribute__ ((aligned(32)));
    Vector2d S_inv_v1  __attribute__ ((aligned(32)));

    Matrix2d S_inv2 __attribute__ ((aligned(32)));
    Vector2d S_inv_v2  __attribute__ ((aligned(32)));

    Matrix2d S_inv3 __attribute__ ((aligned(32)));
    Vector2d S_inv_v3  __attribute__ ((aligned(32)));


    double den[4]  __attribute__ ((aligned(32))), num[4]  __attribute__ ((aligned(32)));

    double vT_S_inv_v[4]  __attribute__ ((aligned(32)));

    __m256d den_v, num_v, vec_vSv, w_vec, zp_dist, zp_angle, z_dist, z_angle, fdiff_dist, fdiff_angle, fdiff1, fdiff2;

    __m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7, ymm8, ymm9, sinv0, sinv1, sinv2, sinv3, s_zeros, s_ones, s_twos, s_threes;

    __m256d determinants, inv_det, neg_inv_det, weights_v;

    __m256i zp_load_mask = _mm256_set_epi64x(2*3*count_zf,2*2*count_zf,2*1*count_zf,0); ///!!!!! Order here
    __m256i sf_mask = _mm256_set_epi64x(4*3*count_zf,4*2*count_zf,4*1*count_zf,0);
    
    for (size_t i = 0; i < NPARTICLES; i+=4) {
        if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
            // std::cout<<"HiFast"<<std::endl;
            //COMPUTE JACOBIANS
            compute_jacobians_fast(particles + i, idf, count_zf, R, zp, Hv, Hf, Sf);
            compute_jacobians_fast(particles + i +1, idf, count_zf, R, zp1, Hv + 1* count_zf, Hf1, Sf1);
            compute_jacobians_fast(particles + i +2, idf, count_zf, R, zp2, Hv + 2* count_zf, Hf2, Sf2);
            compute_jacobians_fast(particles + i +3, idf, count_zf, R, zp3, Hv + 3* count_zf, Hf3, Sf3);
            //END COMPUTE_JACOBIANS

            weights_v = _mm256_load_pd(weights+i);

            for (size_t j = 0; j < count_zf; j++) {
                // Particle 0
                z_dist = _mm256_set1_pd(zf[j][0]);
                z_angle = _mm256_set1_pd(zf[j][1]);

                zp_dist = _mm256_i64gather_pd(*(zp+j), zp_load_mask,8);
                zp_angle = _mm256_i64gather_pd(*(zp+j)+1, zp_load_mask,8);

                fdiff_dist = _mm256_sub_pd(z_dist, zp_dist);
                fdiff_angle = _mm256_sub_pd(z_angle, zp_angle);

                fdiff_angle = simple_pi_to_pi_avx(fdiff_angle);
                
                fdiff1 = _mm256_shuffle_pd(fdiff_dist, fdiff_angle, 0b0000);
                fdiff2 = _mm256_shuffle_pd(fdiff_dist, fdiff_angle, 0b1111);

                _mm256_store2_m128d(feat_diff2[j], feat_diff[j], fdiff1);
                _mm256_store2_m128d(feat_diff3[j], feat_diff1[j], fdiff2);

                s_zeros = _mm256_i64gather_pd(Sf[j], sf_mask,8); //3 1 2 0
                s_ones = _mm256_i64gather_pd(Sf[j]+1, sf_mask,8);
                s_twos = _mm256_i64gather_pd(Sf[j]+2, sf_mask,8);
                s_threes = _mm256_i64gather_pd(Sf[j]+3, sf_mask,8);

                ymm0 = _mm256_mul_pd(s_zeros, s_threes);
                ymm1 = _mm256_mul_pd(s_twos, s_ones);
                
                determinants = _mm256_sub_pd(ymm0, ymm1); //dets0, dets1, dets2, dets3

                inv_det = _mm256_div_pd(_mm256_set1_pd(1.0), determinants);
                neg_inv_det = _mm256_mul_pd(_mm256_set1_pd(-1.0), inv_det);

                s_threes  = _mm256_mul_pd(inv_det, s_threes); //s * A[3]
                s_ones = _mm256_mul_pd(neg_inv_det, s_ones);
                s_twos = _mm256_mul_pd(neg_inv_det, s_twos);
                s_zeros = _mm256_mul_pd(inv_det, s_zeros);

                ymm6 = _mm256_unpacklo_pd(s_threes , s_ones);
                ymm7 = _mm256_unpackhi_pd(s_threes , s_ones);
                ymm8 = _mm256_unpacklo_pd(s_twos, s_zeros);
                ymm9 = _mm256_unpackhi_pd(s_twos, s_zeros);

                _mm256_store_pd(S_inv, _mm256_permute2f128_pd(ymm6, ymm8, 0b00100000));
                _mm256_store_pd(S_inv2,_mm256_permute2f128_pd(ymm6, ymm8, 0b00110001));
                _mm256_store_pd(S_inv1, _mm256_permute2f128_pd(ymm7, ymm9, 0b00100000));
                _mm256_store_pd(S_inv3, _mm256_permute2f128_pd(ymm7, ymm9, 0b00110001));


                // sub(zf[j], zp[j], 2, feat_diff[j]);
                // sub(zf[j], zp1[j], 2, feat_diff1[j]);
                // sub(zf[j], zp2[j], 2, feat_diff2[j]);
                // sub(zf[j], zp3[j], 2, feat_diff3[j]);
                


                // feat_diff[j][1] = pi_to_pi(feat_diff[j][1]);
                // feat_diff1[j][1] = pi_to_pi(feat_diff1[j][1]);
                // feat_diff2[j][1] = pi_to_pi(feat_diff2[j][1]);
                // feat_diff3[j][1] = pi_to_pi(feat_diff3[j][1]);

                // inv_2x2(Sf[j], S_inv);
                // inv_2x2(Sf1[j], S_inv1);
                // inv_2x2(Sf2[j], S_inv2);
                // inv_2x2(Sf3[j], S_inv3);

                mv_2x2(S_inv, feat_diff[j], S_inv_v);
                mv_2x2(S_inv1, feat_diff1[j], S_inv_v1);
                mv_2x2(S_inv2, feat_diff2[j], S_inv_v2);
                mv_2x2(S_inv3, feat_diff3[j], S_inv_v3);

                mul(feat_diff[j], S_inv_v, 1, 2, 1, &vT_S_inv_v[0]); // TODO in linalg   
                mul(feat_diff1[j], S_inv_v1, 1, 2, 1, &vT_S_inv_v[1]); // TODO in linalg
                mul(feat_diff2[j], S_inv_v2, 1, 2, 1, &vT_S_inv_v[2]); // TODO in linalg
                mul(feat_diff3[j], S_inv_v3, 1, 2, 1, &vT_S_inv_v[3]); // TODO in linalg

                KF_cholesky_update(particles[i].xf + 2 * idf[j], particles[i].Pf + 4 * idf[j], 
                                feat_diff[j], R, 
                                Hf[j]);

                KF_cholesky_update(particles[i+1].xf + 2 * idf[j], particles[i+1].Pf + 4 * idf[j], 
                                feat_diff1[j], R, 
                                Hf1[j]);
                KF_cholesky_update(particles[i+2].xf + 2 * idf[j], particles[i+2].Pf + 4 * idf[j], 
                                feat_diff2[j], R, 
                                Hf2[j]);
                KF_cholesky_update(particles[i+3].xf + 2 * idf[j], particles[i+3].Pf + 4 * idf[j], 
                                feat_diff3[j], R, 
                                Hf3[j]);

                // set_xfi(particles+i, particles[i].xf + 2 * idf[j], idf[j]);
                // set_xfi(particles+i+1, particles[i+1].xf + 2 * idf[j], idf[j]);
                // set_xfi(particles+i+2, particles[i+2].xf + 2 * idf[j], idf[j]);
                // set_xfi(particles+i+3, particles[i+3].xf + 2 * idf[j], idf[j]);

                // set_Pfi(particles+i, particles[i].Pf + 4 * idf[j], idf[j]);
                // set_Pfi(particles+i+1, particles[i+1].Pf + 4 * idf[j], idf[j]);
                // set_Pfi(particles+i+2, particles[i+2].Pf + 4 * idf[j], idf[j]);
                // set_Pfi(particles+i+3, particles[i+3].Pf + 4 * idf[j], idf[j]);

                //Weights
                ymm0 = _mm256_sqrt_pd(determinants);
                den_v = _mm256_mul_pd(ymm0, two_pi_vec);

                vec_vSv = _mm256_load_pd(vT_S_inv_v);
                num_v = exp_avx2_pd(_mm256_mul_pd(minus_half, vec_vSv));

                //_mm256_store_pd(num_compare, num_v);
                

                // den[0] = 2 * M_PI * sqrt(determinant_2x2(Sf[j]));
                // num[0] = exp(-0.5 * vT_S_inv_v[0]);
                // w_list[0] *= (double)num[0] / (double)den[0];

                // den[1] = 2 * M_PI * sqrt(determinant_2x2(Sf1[j]));
                // num[1] = exp(-0.5 * vT_S_inv_v[1]);
                // w_list[1] *= (double)num[1] / (double)den[1];

                // den[2] = 2 * M_PI * sqrt(determinant_2x2(Sf2[j]));
                // num[2] = exp(-0.5 * vT_S_inv_v[2]);
                // w_list[2] *= (double)num[2] / (double)den[2];

                // den[3] = 2 * M_PI * sqrt(determinant_2x2(Sf3[j]));
                // num[3] = exp(-0.5 * vT_S_inv_v[3]);
                // w_list[3] *= (double)num[3] / (double)den[3];

                // num[0] = exp(-0.5 * vT_S_inv_v[0]);
                // num[1] = exp(-0.5 * vT_S_inv_v[1]);
                // num[2] = exp(-0.5 * vT_S_inv_v[2]);
                // num[3] = exp(-0.5 * vT_S_inv_v[3]);

                // // for (int i = 0; i<4; i++) {
                // //     std::cout <<num[i]<<", "<<num_compare[i]<<", exp:"<<-0.5 * vT_S_inv_v[i]<<std::endl;
                // //     assert(fabs(num[i]-num_compare[i]) <1.0e-6);
                // // }
                // num_v = _mm256_load_pd(num);
                weights_v = _mm256_mul_pd(weights_v, _mm256_div_pd(num_v, den_v));
                
            }
            _mm256_store_pd(weights+i, weights_v);
                
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature(particles+i, zn, count_zn, R);
            add_feature(particles+i+1, zn, count_zn, R);
            add_feature(particles+i+2, zn, count_zn, R);
            add_feature(particles+i+3, zn, count_zn, R);
        }
    }

    resample_particles(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}

void observe_update_unrolled4x_plain(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    // Seems to be more accurate than base...

    // Compute true data, then add noise
    // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
    // memcpy(ftag_visible, ftag, N_features*sizeof(int));
    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }

    //z is the range and bearing of the observed landmark
    
    get_observations(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
    //print(*z,*Nf_visible,2,std::cout);
    if ( *Nf_visible == 0 ) {
        return;
    }
    
    add_observation_noise(z, *Nf_visible, R, SWITCH_SENSOR_NOISE);

    //Compute (known) data associations
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions
    Vector2d zp[count_zf] __attribute__ ((aligned(32)));
    Matrix23d Hv[4* count_zf] __attribute__ ((aligned(32))); //Unused in whole code
    Matrix2d Hf[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Sf[count_zf] __attribute__ ((aligned(32)));
    Vector2d feat_diff[count_zf] __attribute__ ((aligned(32)));

    Vector2d zp1[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf1[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Sf1[count_zf] __attribute__ ((aligned(32)));
    Vector2d feat_diff1[count_zf] __attribute__ ((aligned(32)));

    Vector2d zp2[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf2[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Sf2[count_zf] __attribute__ ((aligned(32)));
    Vector2d feat_diff2[count_zf] __attribute__ ((aligned(32)));

    Vector2d zp3[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf3[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Sf3[count_zf] __attribute__ ((aligned(32)));
    Vector2d feat_diff3[count_zf] __attribute__ ((aligned(32)));

    //double dx, dy, d2, d, dinv, d2inv, dx_d2inv, dy_d2inv, dx_dinv, dy_dinv;
    //double den, num;
    Matrix2d S_inv;
    Vector2d S_inv_v;

    Matrix2d S_inv1;
    Vector2d S_inv_v1;

    Matrix2d S_inv2;
    Vector2d S_inv_v2;

    Matrix2d S_inv3;
    Vector2d S_inv_v3;


    double den[4], num[4];

    double vT_S_inv_v[4];


    
    for (size_t i = 0; i < NPARTICLES; i+=4) {
        if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
            // std::cout<<"HiFast"<<std::endl;
            //COMPUTE JACOBIANS
            compute_jacobians_fast(particles + i, idf, count_zf, R, zp, Hv, Hf, Sf);
            compute_jacobians_fast(particles + i +1, idf, count_zf, R, zp1, Hv + 1* count_zf, Hf1, Sf1);
            compute_jacobians_fast(particles + i +2, idf, count_zf, R, zp2, Hv + 2* count_zf, Hf2, Sf2);
            compute_jacobians_fast(particles + i +3, idf, count_zf, R, zp3, Hv + 3* count_zf, Hf3, Sf3);
            //END COMPUTE_JACOBIANS
            double w_list[4] = {1.0,1.0,1.0,1.0};

            for (size_t j = 0; j < count_zf; j++) {
                // Particle 0
                sub(zf[j], zp[j], 2, feat_diff[j]);
                sub(zf[j], zp1[j], 2, feat_diff1[j]);
                sub(zf[j], zp2[j], 2, feat_diff2[j]);
                sub(zf[j], zp3[j], 2, feat_diff3[j]);

                feat_diff[j][1] = pi_to_pi(feat_diff[j][1]);
                feat_diff1[j][1] = pi_to_pi(feat_diff1[j][1]);
                feat_diff2[j][1] = pi_to_pi(feat_diff2[j][1]);
                feat_diff3[j][1] = pi_to_pi(feat_diff3[j][1]);


                inv_2x2(Sf[j], S_inv);
                inv_2x2(Sf1[j], S_inv1);
                inv_2x2(Sf2[j], S_inv2);
                inv_2x2(Sf3[j], S_inv3);

                mv_2x2(S_inv, feat_diff[j], S_inv_v);
                mv_2x2(S_inv, feat_diff1[j], S_inv_v1);
                mv_2x2(S_inv2, feat_diff2[j], S_inv_v2);
                mv_2x2(S_inv3, feat_diff3[j], S_inv_v3);

                mul(feat_diff[j], S_inv_v, 1, 2, 1, &vT_S_inv_v[0]); // TODO in linalg   
                mul(feat_diff1[j], S_inv_v1, 1, 2, 1, &vT_S_inv_v[1]); // TODO in linalg
                mul(feat_diff2[j], S_inv_v2, 1, 2, 1, &vT_S_inv_v[2]); // TODO in linalg
                mul(feat_diff3[j], S_inv_v3, 1, 2, 1, &vT_S_inv_v[3]); // TODO in linalg

                KF_cholesky_update(particles[i].xf + 2 * idf[j], particles[i].Pf + 4 * idf[j], 
                                feat_diff[j], R, 
                                Hf[j]);
                KF_cholesky_update(particles[i+1].xf + 2 * idf[j], particles[i+1].Pf + 4 * idf[j], 
                                feat_diff1[j], R, 
                                Hf1[j]);
                KF_cholesky_update(particles[i+2].xf + 2 * idf[j], particles[i+2].Pf + 4 * idf[j], 
                                feat_diff2[j], R, 
                                Hf2[j]);
                KF_cholesky_update(particles[i+3].xf + 2 * idf[j], particles[i+3].Pf + 4 * idf[j], 
                                feat_diff3[j], R, 
                                Hf3[j]);

                set_xfi(particles+i, particles[i].xf + 2 * idf[j], idf[j]);
                set_xfi(particles+i+1, particles[i+1].xf + 2 * idf[j], idf[j]);
                set_xfi(particles+i+2, particles[i+2].xf + 2 * idf[j], idf[j]);
                set_xfi(particles+i+3, particles[i+3].xf + 2 * idf[j], idf[j]);

                set_Pfi(particles+i, particles[i].Pf + 4 * idf[j], idf[j]);
                set_Pfi(particles+i+1, particles[i+1].Pf + 4 * idf[j], idf[j]);
                set_Pfi(particles+i+2, particles[i+2].Pf + 4 * idf[j], idf[j]);
                set_Pfi(particles+i+3, particles[i+3].Pf + 4 * idf[j], idf[j]);

                //Weights
                den[0] = 2 * M_PI * sqrt(determinant_2x2(Sf[j]));
                num[0] = exp(-0.5 * vT_S_inv_v[0]);
                w_list[0] *= (double)num[0] / (double)den[0];

                den[1] = 2 * M_PI * sqrt(determinant_2x2(Sf1[j]));
                num[1] = exp(-0.5 * vT_S_inv_v[1]);
                w_list[1] *= (double)num[1] / (double)den[1];

                den[2] = 2 * M_PI * sqrt(determinant_2x2(Sf2[j]));
                num[2] = exp(-0.5 * vT_S_inv_v[2]);
                w_list[2] *= (double)num[2] / (double)den[2];

                den[3] = 2 * M_PI * sqrt(determinant_2x2(Sf3[j]));
                num[3] = exp(-0.5 * vT_S_inv_v[3]);
                w_list[3] *= (double)num[3] / (double)den[3];

                
            }
            weights[i]*=w_list[0];
            weights[i+1]*=w_list[1];
            weights[i+2]*=w_list[2];
            weights[i+3]*=w_list[3];
                
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature(particles+i, zn, count_zn, R);
            add_feature(particles+i+1, zn, count_zn, R);
            add_feature(particles+i+2, zn, count_zn, R);
            add_feature(particles+i+3, zn, count_zn, R);
        }
    }

    resample_particles(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}

void observe_update_unrolled4x(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    // Exactly like base.

    // Compute true data, then add noise
    // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
    // memcpy(ftag_visible, ftag, N_features*sizeof(int));
    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }

    //z is the range and bearing of the observed landmark
    
    get_observations(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
    //print(*z,*Nf_visible,2,std::cout);
    if ( *Nf_visible == 0 ) {
        return;
    }
    
    add_observation_noise(z, *Nf_visible, R, SWITCH_SENSOR_NOISE);

    //Compute (known) data associations
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions
    Vector2d zp[4* count_zf] __attribute__ ((aligned(32)));
    Matrix23d Hv[4* count_zf] __attribute__ ((aligned(32))); //Unused in whole code
    Matrix2d Hf[4* count_zf] __attribute__ ((aligned(32)));
    Matrix2d Sf[4 * count_zf] __attribute__ ((aligned(32)));

    Vector2d feat_diff[4 * count_zf] __attribute__ ((aligned(32)));  // difference btw feature prediciton and
                                        // measurement (used to update mean)
    // perform update

    //double dx, dy, d2, d, dinv, d2inv, dx_d2inv, dy_d2inv, dx_dinv, dy_dinv;
    double den, num;
    Matrix2d S_inv;
    Vector2d S_inv_v;
    double vT_S_inv_v;
    for (size_t i = 0; i < NPARTICLES; i+=4) {
        if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
            // std::cout<<"HiFast"<<std::endl;
            //COMPUTE JACOBIANS
            compute_jacobians_fast(particles + i, idf, count_zf, R, zp, Hv, Hf, Sf);
            compute_jacobians_fast(particles + i +1, idf, count_zf, R, zp + 1* count_zf, Hv + 1* count_zf, Hf + 1* count_zf, Sf + 1* count_zf);
            compute_jacobians_fast(particles + i +2, idf, count_zf, R, zp + 2* count_zf, Hv + 2* count_zf, Hf + 2* count_zf, Sf + 2* count_zf);
            compute_jacobians_fast(particles + i +3, idf, count_zf, R, zp + 3* count_zf, Hv + 3* count_zf, Hf + 3* count_zf, Sf + 3* count_zf);
            //END COMPUTE_JACOBIANS
            for (int k = 0; k<4; k++) {
                // std::cerr<<std::endl<<"k="<<k<<", i="<<i<<std::endl;
                double w = 1.0;
                for (size_t j = 0; j < count_zf; j++) {

                    sub(zf[j], zp[j + k* count_zf], 2, feat_diff[j + k* count_zf]);
                    feat_diff[j + k* count_zf][1] = pi_to_pi(feat_diff[j + k* count_zf][1]);
                    inv_2x2(Sf[j + k* count_zf], S_inv);
                    mv_2x2(S_inv, feat_diff[j + k* count_zf], S_inv_v);
                    mul(feat_diff[j + k* count_zf], S_inv_v, 1, 2, 1, &vT_S_inv_v); // TODO in linalg
                    den = 2 * M_PI * sqrt(determinant_2x2(Sf[j + k* count_zf]));
                    num = exp(-0.5 * vT_S_inv_v);
                    w *= (double)num / (double)den;

                    
                    KF_cholesky_update_base(particles[i+k].xf + 2 * idf[j], particles[i + k].Pf + 4 * idf[j], 
                                    feat_diff[j + k* count_zf], R, Hf[j + k* count_zf]);
                    set_xfi(particles+i+k, particles[i+k].xf + 2 * idf[j], idf[j]);
                    set_Pfi(particles+i+k, particles[i+k].Pf + 4 * idf[j], idf[j]);
                }
                weights[i+k]*=w;
            }
                
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature(particles+i, zn, count_zn, R);
            add_feature(particles+i+1, zn, count_zn, R);
            add_feature(particles+i+2, zn, count_zn, R);
            add_feature(particles+i+3, zn, count_zn, R);
        }
    }

    resample_particles(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}

void observe_update_simplified(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    // Compute true data, then add noise
    // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
    // memcpy(ftag_visible, ftag, N_features*sizeof(int));
    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }

    //z is the range and bearing of the observed landmark
    
    get_observations(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
    //print(*z,*Nf_visible,2,std::cout);
    if ( *Nf_visible == 0 ) {
        return;
    }
    
    add_observation_noise(z, *Nf_visible, R, SWITCH_SENSOR_NOISE);

    //Compute (known) data associations
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions
    Vector2d zp[count_zf] __attribute__ ((aligned(32)));
    Matrix23d Hv[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Sf[count_zf] __attribute__ ((aligned(32)));

    Vector2d v[count_zf]  __attribute__ ((aligned(32)));
    Vector2d feat_diff[count_zf];  // difference btw feature prediciton and
                                        // measurement (used to update mean)
    // perform update
    //double dx, dy, d2, d, dinv, d2inv, dx_d2inv, dy_d2inv, dx_dinv, dy_dinv;
    double den, num;
    Matrix2d S, ST, S_inv;
    Vector2d S_inv_v;
    double vT_S_inv_v;
    for (size_t i = 0; i < NPARTICLES; i++) {
        if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
            // std::cout<<"HiFast"<<std::endl;
            //COMPUTE JACOBIANS
            compute_jacobians_fast(particles + i, idf, count_zf, R, zp, Hv, Hf, Sf);
            //END COMPUTE_JACOBIANS
            double w = 1.0;
            for (size_t j = 0; j < count_zf; j++) {
                sub(zf[j], zp[j], 2, feat_diff[j]);
                feat_diff[j][1] = pi_to_pi(feat_diff[j][1]); //=v[j]
                
                copy(Sf[j],4,S);
                inv_2x2(Sf[j], S_inv);
                
                mv_2x2(S_inv, feat_diff[j], S_inv_v);
                mul(feat_diff[j], S_inv_v, 1, 2, 1, &vT_S_inv_v); // TODO in linalg

                den = 2 * M_PI * sqrt(determinant_2x2(Sf[j]));
                num = exp(-0.5 * vT_S_inv_v);
                w *= (double)num / (double)den;

                
                KF_cholesky_update_base(particles[i].xf + 2 * idf[j], particles[i].Pf + 4 * idf[j], 
                                feat_diff[j], R, 
                                Hf[j]);
                //set_xfi(particles+i, particles[i].xf + 2 * idf[j], idf[j]);
                //set_Pfi(particles+i, particles[i].Pf + 4 * idf[j], idf[j]);
                
            }
            weights[i]*=w;
            
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature(particles+i, zn, count_zn, R);
        }
    }

    resample_particles(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}

void observe_update_inplace(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    // Compute true data, then add noise
    // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
    // memcpy(ftag_visible, ftag, N_features*sizeof(int));
    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }

    //z is the range and bearing of the observed landmark
    
    get_observations(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
    //print(*z,*Nf_visible,2,std::cout);
    if ( *Nf_visible == 0 ) {
        return;
    }
    
    add_observation_noise(z, *Nf_visible, R, SWITCH_SENSOR_NOISE);

    //Compute (known) data associations
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions

    // perform update
    for (size_t i = 0; i < NPARTICLES; i++) {
        if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
            Vector2d zp[count_zf] __attribute__ ((aligned(32)));
            Matrix23d Hv[count_zf] __attribute__ ((aligned(32)));
            Matrix2d Hf[count_zf] __attribute__ ((aligned(32)));
            Matrix2d Sf[count_zf] __attribute__ ((aligned(32)));
            compute_jacobians_fast(&particles[i], idf, count_zf, R, zp, Hv, Hf, Sf);
            
            Vector2d v[count_zf];
            for (size_t j = 0; j < count_zf; j++) {
                Vector2d v_j;
                sub(zf[j], zp[j], 2, v_j);  // v_j = z[j] - zp[j]
                v_j[1] = pi_to_pi_base(v_j[1]);
                copy(v_j, 2, v[j]);  // v[j] = v_j
            }

            double w0 = 1.0;

            double den, num;
            // this can probably be done alot faster without this loop.....
            
            for (size_t j = 0; j <count_zf; j++) {
                Matrix2d S, ST, S_inv;
                Vector2d S_inv_v;
                double vT_S_inv_v;
                // Eq. 61 in Thrun03g
                
                copy(Sf[j], 4, S);
                transpose(S, 2, 2, ST);
                inv_2x2(S, S_inv);

                
                mv_2x2(S_inv, v[j], S_inv_v);
                mul(v[j], S_inv_v, 1, 2, 1, &vT_S_inv_v); // TODO in linalg

                print(v[j],2,1);
                print(S_inv,4,1);
                assert(num);
                den = 2 * M_PI * sqrt(determinant_2x2(S));
                num = exp(-0.5 * vT_S_inv_v);
                w0 *= (double)num / (double)den;
            }
            double w = w0;
            
            //double w = compute_weight_active(&particles[i], zf, count_zf, idf, R, zp, Hv, Hf, Sf);
            w*=weights[i];
            weights[i] = w;
            feature_update_active(&particles[i], zf, idf, count_zf, R, zp, Hv, Hf, Sf);
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature(&particles[i], zn, count_zn, R);
        }
    }

    resample_particles(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}



