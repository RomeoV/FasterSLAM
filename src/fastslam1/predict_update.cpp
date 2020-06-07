#include "predict_update.h"

#include "particle.h"

#include "compute_steering.h"
#include "predict_true.h"
#include "add_control_noise.h"
#include "predict.h"
#include "configfile.h"
#include "pi_to_pi.h"
#include "trigonometry.h"
#include <math.h>
#include "tscheb_sine.h"
#include "linalg.h"




void predict_update(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
#ifdef __AVX2__
    predict_update_fast(wp, N_waypoints, V, Q, dt, N, xtrue, iwp, G, particles);                   
#else
#warning "Using predict_update_base because AVX2 is not supported!"
    predict_update_base(wp, N_waypoints, V, Q, dt, N, xtrue, iwp, G, particles);                   
#endif
}

//! --------------- !//
//! --- VP BASE --- !//
//! --------------- !//

void predict_VP(Vector3d state, double V_, double G_, double *Q, double WB, double dt, bool add_control_noise) {
    V_ = V_/(1.0-tan(G_) *0.76/2.83);
    if (add_control_noise) {
        double VG[2] = {V_,G_};
        double VnGn[2];
        multivariate_gauss_base(VG,Q,VnGn);
        V_=VnGn[0];
        G_=VnGn[1];
    }

    double a=3.78;
    double b=0.5;
    
    state[0]+= V_*dt*cos(state[2]) - V_ /WB *tan(G_) * dt* (a*sin(state[2]) + b * cos(state[2]));
    state[1]+= V_*dt*sin(state[2]) + V_ /WB *tan(G_) * dt* (a*cos(state[2]) - b * sin(state[2]));
    state[2] = pi_to_pi_base(state[2] + V_*dt*tan(G_)/WB);
}

void predict_update_VP_base(double* controls, size_t N_controls, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    predict_VP(xtrue,V,*G, Q, WHEELBASE, dt, false);

    double VnGn[2];
    add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    for (size_t i = 0; i < N; i++) {
        predict_VP(particles[i].xv, VnGn[0], VnGn[1], Q, WHEELBASE, dt,SWITCH_PREDICT_NOISE);
    }
}

//! --------------- !//
//! -- VP ACTIVE -- !//
//! --------------- !//

void predict_VP_active(Vector3d state, double V_, double G_, double *S, double WB, double dt, bool add_control_noise) {

    V_ = V_ / ( 1.0 - tan(G_)*0.76/2.83 );
    
    if (add_control_noise) {
        double X[2];
        fill_rand(X, 2, -1.0, 1.0);
        
        double VnGn[2] = {V_, G_};    
        mvadd_2x2(S, X, VnGn);

        V_ = VnGn[0];
        G_ = VnGn[1];
    }

    double a = 3.78;
    double b = 0.5;
    const double V_dt = V_*dt;
    const double alpha = V_dt * tan(G_) / WB;
    const double cos2 = cos(state[2]);
    const double sin2 = sin(state[2]);
    const double alpha_a = alpha*a;
    const double alpha_b = alpha*b;
    const double beta = V_dt - alpha_b;

    state[0] += beta*cos2 - alpha_a*sin2;
    state[1] += beta*sin2 + alpha_a*cos2;
    state[2] = pi_to_pi_active(state[2] + alpha);
}

void predict_VP_unrolledx4_active(Vector3d state0, Vector3d state1, Vector3d state2, Vector3d state3, 
                                    double V_, double G_, double *S, double WB, double dt, bool add_control_noise) {
    double V0_ = V_;
    double V1_ = V_;
    double V2_ = V_;
    double V3_ = V_;

    double G0_ = G_;
    double G1_ = G_;
    double G2_ = G_;
    double G3_ = G_;

    if (add_control_noise) { 
        double VnGn01[4] = {V_, G_, V_, G_};
        double VnGn23[4] = {V_, G_, V_, G_};
         
        double X01[4], X23[4];
        fill_rand(X01, 4, -1.0, 1.0);
        fill_rand(X23, 4, -1.0, 1.0);
        
        double ST[4] = { };
        transpose_2x2(S, ST);
        mmadd_2x2(X01, ST, VnGn01); // S transpose 
        mmadd_2x2(X23, ST, VnGn23); // S transpose
 
        V0_ = VnGn01[0];
        G0_ = VnGn01[1];
        V1_ = VnGn01[2];
        G1_ = VnGn01[3];
        
        V2_ = VnGn23[0];
        G2_ = VnGn23[1];
        V3_ = VnGn23[2];
        G3_ = VnGn23[3];
    }

    const double a = 3.78;
    const double b = 0.5;

    const double V0_dt = V0_*dt;
    const double V1_dt = V1_*dt;
    const double V2_dt = V2_*dt;
    const double V3_dt = V3_*dt;

    const double alpha0 = V0_dt * tan( G0_ ) / WB;
    const double alpha1 = V1_dt * tan( G1_ ) / WB;
    const double alpha2 = V2_dt * tan( G2_ ) / WB;
    const double alpha3 = V3_dt * tan( G3_ ) / WB;
    
    const double cos02 = cos( state0[2] );
    const double cos12 = cos( state1[2] );
    const double cos22 = cos( state2[2] );
    const double cos32 = cos( state3[2] );

    const double sin02 = sin( state0[2] );
    const double sin12 = sin( state1[2] );
    const double sin22 = sin( state2[2] );
    const double sin32 = sin( state3[2] );

    const double alpha0_a = alpha0*a;
    const double alpha1_a = alpha1*a;
    const double alpha2_a = alpha2*a;
    const double alpha3_a = alpha3*a;
    
    const double alpha0_b = alpha0*b;
    const double alpha1_b = alpha1*b;
    const double alpha2_b = alpha2*b;
    const double alpha3_b = alpha3*b;
    
    const double beta0 = V0_dt - alpha0_b;
    const double beta1 = V1_dt - alpha1_b;
    const double beta2 = V2_dt - alpha2_b;
    const double beta3 = V3_dt - alpha3_b;

    state0[0] += beta0*cos02 - alpha0_a*sin02;
    state1[0] += beta1*cos12 - alpha1_a*sin12;
    state2[0] += beta2*cos22 - alpha2_a*sin22;
    state3[0] += beta3*cos32 - alpha3_a*sin32;
    
    state0[1] += beta0*sin02 + alpha0_a*cos02;
    state1[1] += beta1*sin12 + alpha1_a*cos12;
    state2[1] += beta2*sin22 + alpha2_a*cos22;
    state3[1] += beta3*sin32 + alpha3_a*cos32;
    
    state0[2] = pi_to_pi_active(state0[2] + alpha0);
    state1[2] = pi_to_pi_active(state1[2] + alpha1);
    state2[2] = pi_to_pi_active(state2[2] + alpha2);
    state3[2] = pi_to_pi_active(state3[2] + alpha3);
}

void predict_VP_unrolledx4_active_avx(Vector3d state0,
                                      Vector3d state1,
                                      Vector3d state2,
                                      Vector3d state3,
                                      double V_,
                                      double G_,
                                      double *S,
                                      double WB,
                                      double dt,
                                      bool add_control_noise) {

    double const a = 3.78;
    double const b = 0.5;
    __m256d const as = _mm256_set1_pd( a );
    __m256d const bs = _mm256_set1_pd( b );
    __m256d const dts = _mm256_set1_pd( dt ); 
    __m256d const Ss = _mm256_load_pd( S );
    __m256d const invWB = _mm256_set1_pd( 1.0 / WB );
    __m256d Vs = _mm256_set1_pd( V_ );
    __m256d Gs = _mm256_set1_pd( G_ );
    __m256d state_0 = _mm256_set_pd( state3[0], state2[0], state1[0], state0[0] );
    __m256d state_1 = _mm256_set_pd( state3[1], state2[1], state1[1], state0[1] );
    __m256d state_2 = _mm256_set_pd( state3[2], state2[2], state1[2], state0[2] );

    if (add_control_noise) { 
        __m256d const VnGn = _mm256_set_pd(G_, V_, G_, V_);

        double X02[4] __attribute__ ((aligned(32)));
        double X13[4] __attribute__ ((aligned(32)));
        fill_rand(X02+0, 2, -1.0, 1.0);
        fill_rand(X13+0, 2, -1.0, 1.0);
        fill_rand(X02+2, 2, -1.0, 1.0);
        fill_rand(X13+2, 2, -1.0, 1.0);
        __m256d const x02 = _mm256_load_pd( X02 );
        __m256d const x13 = _mm256_load_pd( X13 );
//        __m256d const x02 = fill_rand_avx(-1.0, 1.0);
//        __m256d const x13 = fill_rand_avx(-1.0, 1.0);
        
        __m256d const VnGn02 = _mmTadd_2x2_avx_v2(x02, Ss, VnGn);
        __m256d const VnGn13 = _mmTadd_2x2_avx_v2(x13, Ss, VnGn);
 
        Vs = _mm256_shuffle_pd(VnGn02, VnGn13, 0b0000);
        Gs = _mm256_shuffle_pd(VnGn02, VnGn13, 0b1111);
    }

    __m256d const Vsdt  = _mm256_mul_pd(Vs, dts);
    __m256d const sinGs = tscheb_sin_avx( Gs );
    __m256d const cosGs = tscheb_cos_avx( Gs );
    __m256d const tanGs = _mm256_div_pd( sinGs, cosGs );
    __m256d const alpha = _mm256_mul_pd( _mm256_mul_pd( Vsdt, tanGs ), invWB );
  
    __m256d const sin_state_2 = tscheb_sin_avx( state_2 );
    __m256d const cos_state_2 = tscheb_cos_avx( state_2 );

    __m256d const alpha_a = _mm256_mul_pd(alpha, as);
    __m256d const alpha_b = _mm256_mul_pd(alpha, bs);
    __m256d const beta    = _mm256_sub_pd(Vsdt, alpha_b);
    
    state_0 = _mm256_fmadd_pd(beta,    cos_state_2, state_0);
    state_0 = _mm256_sub_pd( state_0, _mm256_mul_pd(alpha_a, sin_state_2));
 
    state_1 = _mm256_fmadd_pd(beta,    sin_state_2, state_1);
    state_1 = _mm256_fmadd_pd(alpha_a, cos_state_2, state_1);
    
    state_2 = simple_pi_to_pi_avx(_mm256_add_pd(state_2, alpha));

    __m256d state_3 = _mm256_setzero_pd(); // dummy 
    
    __m256d t0, t1, t2, t3;
    register_transpose(state_0, state_1, state_2, state_3, &t0, &t1, &t2, &t3);

    __m256i mask = _mm256_set_epi64x(+1, -1, -1, -1); // store only the first 3
    _mm256_maskstore_pd(state0, mask, t0);
    _mm256_maskstore_pd(state1, mask, t1);
    _mm256_maskstore_pd(state2, mask, t2);
    _mm256_maskstore_pd(state3, mask, t3);
}

// __m256d tscheb_sin_avx(__m256d alphas)
// __m256d tscheb_cos_avx(__m256d alphas)
void predict_update_VP_active(double* controls, size_t N_controls, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    
    double S[4] __attribute__((aligned(32))) = {};
    llt_2x2(Q, S);
    
    predict_VP_active(xtrue, V, *G, S, WHEELBASE, dt, false); // Pass Cholesky factor S instead of Q

    double VnGn[2];
    add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn);
 
    double const a = 3.78;
    double const b = 0.5;
    __m256d const as = _mm256_set1_pd( a );
    __m256d const bs = _mm256_set1_pd( b );
    
    __m256d const dts = _mm256_set1_pd( dt ); 
    
    __m256d const Ss = _mm256_load_pd( S );
    
    __m256d const invWB = _mm256_set1_pd( 1.0 / WHEELBASE);
    
    VnGn[0] = VnGn[0] / ( 1.0 - tan( VnGn[1] )*0.76/2.83 ); // predict_VP_unrolledx4 takes as input the transformed VnGn[0] ( reuse )
    __m256d Vs = _mm256_set1_pd( VnGn[0] );
    __m256d Gs = _mm256_set1_pd( VnGn[1] );
    __m256d const VnGnv = _mm256_set_pd(VnGn[1], VnGn[0], VnGn[1], VnGn[0]);
        
    __m256d state_3 = _mm256_setzero_pd(); // dummy 
    __m256i mask = _mm256_set_epi64x(+1, -1, -1, -1); // store only the first 3
    
    for (size_t i = 0; i < N; i+=4) {
        //predict_VP_unrolledx4_active_avx(particles[i+0].xv,
        //                                 particles[i+1].xv,
        //                                 particles[i+2].xv,
        //                                 particles[i+3].xv,
        //                                 VnGn[0], VnGn[1], S, // Pass Cholesky factor S instead of Q for reuse
        //                                 WHEELBASE, dt, SWITCH_PREDICT_NOISE);
        // Load transposed states
 
        //__m256i mask = _mm256_set_epi64x(+1, -1, -1, -1); // store only the first 3
        //__m256d p0 = _mm256_maskload_pd(particles[i+0].xv, mask);
        //__m256d p1 = _mm256_maskload_pd(particles[i+1].xv, mask);
        //__m256d p2 = _mm256_maskload_pd(particles[i+2].xv, mask);
        //__m256d p3 = _mm256_maskload_pd(particles[i+3].xv, mask);
        //__m256d state_0, state_1, state_2, state_3;
        //register_transpose(p0, p1, p2, p3, &state_0, &state_1, &state_2, &state_3);

        __m256d state_0 = _mm256_set_pd( particles[i+3].xv[0], particles[i+2].xv[0], particles[i+1].xv[0], particles[i+0].xv[0] );
        __m256d state_1 = _mm256_set_pd( particles[i+3].xv[1], particles[i+2].xv[1], particles[i+1].xv[1], particles[i+0].xv[1] );
        __m256d state_2 = _mm256_set_pd( particles[i+3].xv[2], particles[i+2].xv[2], particles[i+1].xv[2], particles[i+0].xv[2] );

        if (add_control_noise) { 

            //double X02[4] __attribute__ ((aligned(32)));
            //double X13[4] __attribute__ ((aligned(32)));
            //fill_rand(X02+0, 2, -1.0, 1.0);
            //fill_rand(X13+0, 2, -1.0, 1.0);
            //fill_rand(X02+2, 2, -1.0, 1.0);
            //fill_rand(X13+2, 2, -1.0, 1.0);
            //__m256d const x02 = _mm256_load_pd( X02 );
            //__m256d const x13 = _mm256_load_pd( X13 );
            __m256d const x02 = fill_rand_avx(-1.0, 1.0);
            __m256d const x13 = fill_rand_avx(-1.0, 1.0);

            __m256d const VnGn02 = _mmTadd_2x2_avx_v2(x02, Ss, VnGnv);
            __m256d const VnGn13 = _mmTadd_2x2_avx_v2(x13, Ss, VnGnv);

            Vs = _mm256_shuffle_pd(VnGn02, VnGn13, 0b0000);
            Gs = _mm256_shuffle_pd(VnGn02, VnGn13, 0b1111);
        }

        __m256d const Vsdt  = _mm256_mul_pd(Vs, dts);
        __m256d const sinGs = tscheb_sin_avx( Gs );
        __m256d const cosGs = tscheb_cos_avx( Gs );
        __m256d const tanGs = _mm256_div_pd( sinGs, cosGs );
        __m256d const alpha = _mm256_mul_pd( _mm256_mul_pd( Vsdt, tanGs ), invWB );

        __m256d const sin_state_2 = tscheb_sin_avx( state_2 );
        __m256d const cos_state_2 = tscheb_cos_avx( state_2 );

        __m256d const alpha_a = _mm256_mul_pd(alpha, as);
        __m256d const alpha_b = _mm256_mul_pd(alpha, bs);
        __m256d const beta    = _mm256_sub_pd(Vsdt, alpha_b);

        state_0 = _mm256_fmadd_pd(beta,    cos_state_2, state_0);
        state_0 = _mm256_sub_pd( state_0, _mm256_mul_pd(alpha_a, sin_state_2));

        state_1 = _mm256_fmadd_pd(beta,    sin_state_2, state_1);
        state_1 = _mm256_fmadd_pd(alpha_a, cos_state_2, state_1);

        state_2 = simple_pi_to_pi_avx(_mm256_add_pd(state_2, alpha));

        //__m256d state_3 = _mm256_setzero_pd(); // dummy 

        __m256d t0, t1, t2, t3;
        register_transpose(state_0, state_1, state_2, state_3, &t0, &t1, &t2, &t3);

        //__m256i mask = _mm256_set_epi64x(+1, -1, -1, -1); // store only the first 3
        _mm256_maskstore_pd(particles[i+0].xv, mask, t0);
        _mm256_maskstore_pd(particles[i+1].xv, mask, t1);
        _mm256_maskstore_pd(particles[i+2].xv, mask, t2);
        _mm256_maskstore_pd(particles[i+3].xv, mask, t3);
    }
}

/*****************************************************************************
 * PERFORMANCE STATUS (N=NPARTICLES)
 * Work, best: 20 (pred_true) + 22(com_steer.) + N*20 (predict) = 42 + N*20 flops
 * Work, worst: 29 (pred_true) +17 (con_noise) + 37(com_steer.) + N*46 (predict) = 83 + N*46 flops
 * 
 * 
 * #work best, det: 1 atan2 + 2 pow_2 + 3*N+3 sin/cos + 9*N+9+3 mults + 3*N+3+8 adds + 1*N+1+1 neg + 4*N+4+7 fl-comp + 1*N+1 div  
 * #work, worst, det: 1 atan2 + 2 pow_2 + 3*N+3+2 neg + 17*N+12+9 mults + 4*N+2+1 div + 1*N+1+1 floor + 4*N+4+11 fl-comp + 3*N+3 sin + 12*N+5+10 adds + 2*N sqrt 
 * #Work best, detailed: 
 * Memory moved: TBD
 * Cycles: 48000 with N = 100
 * Performance: 0.04
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void predict_update_base(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    compute_steering_base(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
            if ( *iwp == -1 && NUMBER_LOOPS > 1 ) {
                *iwp = 0;
                NUMBER_LOOPS--;
            }
    predict_true_base(V, *G, WHEELBASE, dt, xtrue);

    // add process noise
    double VnGn[2];
    add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    // Predict step	
    for (size_t i = 0; i < N; i++) {
        predict_base(&particles[i], VnGn[0], VnGn[1], Q, WHEELBASE, dt);
    }
}

double predict_update_base_flops(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    int iwp_reset = *iwp;
    double G_reset = *G;
    double flop_count = compute_steering_base_flops(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
    
    flop_count += predict_true_base_flops(V, *G, WHEELBASE, dt, xtrue);

    double VnGn[2] = {0,0};
    flop_count +=add_control_noise_base_flops(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn);
    //add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn);

    //flop_count += add_control_noise_base_flops(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    for (size_t i = 0; i < N; i++) {
        flop_count += predict_base_flops(&particles[i], VnGn[0], VnGn[1], Q, WHEELBASE, dt);
    }
    *iwp = iwp_reset;
    *G = G_reset;
    return flop_count;
}

double predict_update_base_memory(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    int iwp_reset = *iwp;
    double G_reset = *G;
    double memory_moved = compute_steering_base_memory(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
    memory_moved += predict_true_base_memory(V, *G, WHEELBASE, dt, xtrue);
    double VnGn[2] = {0,0};
    //add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn);
    //memory_moved += add_control_noise_base_memory(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    for (size_t i = 0; i < N; i++) {
        memory_moved += predict_base_memory(&particles[i], VnGn[0], VnGn[1], Q, WHEELBASE, dt);
    }
    return memory_moved;
}




void predict_update_active(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    compute_steering(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
            if ( *iwp == -1 && NUMBER_LOOPS > 1 ) {
                *iwp = 0;
                NUMBER_LOOPS--;
            }
    predict_true(V, *G, WHEELBASE, dt, xtrue);

    // add process noise
    double VnGn[2];
    add_control_noise(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    // Predict step	
    for (size_t i = 0; i < N; i++) {
        predict(&particles[i], VnGn[0], VnGn[1], Q, WHEELBASE, dt);
    }
}


#ifdef __AVX2__
void predict_update_fast_plain(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
                        
    compute_steering_base(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
            if ( *iwp == -1 && NUMBER_LOOPS > 1 ) {
                *iwp = 0;
                NUMBER_LOOPS--;
            }
    predict_true_base(V, *G, WHEELBASE, dt, xtrue);

    // add process noise
    double VnGn[2] __attribute__ ((aligned(32)));

    add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    // Predict step
    double Vn, Gn, Vn1, Gn1, Vn2, Gn2, Vn3, Gn3;
    __m256d Vns, VndtWB, Vndt;
    __m256d Gns, Gn_theta, Gnsin;

    __m256d angles;

    __m256d dtv = _mm256_set1_pd(dt);
    __m256d WBv = _mm256_set1_pd(WHEELBASE);
    __m256d dtWBv = _mm256_div_pd(dtv, WBv);

    __m256i load_mask = _mm256_set_epi64x(9,6,3,0);

    __m256d thetas, xs, ys, xv1, xv2, xv3, xv, Gn_theta_sin, Gn_theta_cos;
    Vns = _mm256_set1_pd(VnGn[0]);
    Gns = _mm256_set1_pd(VnGn[1]);
    double angle_buffer[4] __attribute__ ((aligned(32)));

    double S[4]  __attribute__ ((aligned(32)));

    llt_2x2(Q,S);
    
    __m256d SMat = _mm256_load_pd(S);
    __m256d VG = _mm256_set_pd(VnGn[1], VnGn[0], VnGn[1], VnGn[0]);


    for (size_t i = 0; i < N; i+=4) {
        if (SWITCH_PREDICT_NOISE == 1) {
            //Vector2d noise;
            //multivariate_gauss_base(VnGn,Q,noise);	
            __m256d rand_vec1 = fill_rand_avx(-1.0,1.0);
            __m256d rand_vec2 = fill_rand_avx(-1.0,1.0);

            __m256d mul_gauss1 = _mmTadd_2x2_avx_v2(rand_vec1, SMat, VG);
            __m256d mul_gauss2 = _mmTadd_2x2_avx_v2(rand_vec2, SMat, VG);

            // S*rand + VnGn

            Vns = _mm256_shuffle_pd(mul_gauss1, mul_gauss2,0b0000);
            Gns = _mm256_shuffle_pd(mul_gauss1, mul_gauss2,0b1111);
        }
        xs = _mm256_i64gather_pd(particles[i].xv, load_mask, 8);
        ys = _mm256_i64gather_pd(particles[i].xv + 1, load_mask, 8);

        thetas = _mm256_i64gather_pd(particles[i].xv +2, load_mask, 8);

        Vndt = _mm256_mul_pd(Vns, dtv);
        VndtWB = _mm256_mul_pd(Vns, dtWBv);
        Gn_theta = _mm256_add_pd(Gns, thetas);
        Gn_theta_cos = tscheb_cos_avx(Gn_theta);
        Gn_theta_sin = tscheb_sin_avx(Gn_theta);
        Gnsin = tscheb_sin_avx(Gns);


        xs = _mm256_fmadd_pd(Vndt, Gn_theta_cos, xs);
        ys = _mm256_fmadd_pd(Vndt, Gn_theta_sin, ys);
        angles = _mm256_fmadd_pd(VndtWB, Gnsin, thetas);
        
        xv1 = _mm256_shuffle_pd(xs,ys,0b0000); // x1, y1, x3, y3
        xv2 = _mm256_shuffle_pd(xs,ys,0b1111); // x2, y2, x4, y4

        _mm256_store2_m128d(particles[i+2].xv, particles[i].xv, xv1);
        _mm256_store2_m128d(particles[i+3].xv, particles[i+1].xv, xv2);
        
        angles = simple_pi_to_pi_avx(angles);
        _mm256_store_pd(angle_buffer, angles);
        //double xv2 = particles[i].xv[2];
        //particles[i].xv[0] += Vn*dt*tscheb_cos(Gn + xv2);
        //particles[i].xv[1] += Vn*dt*tscheb_sin(Gn + xv2); 

        for (int j = 0; j<4; j++) {
            particles[i+j].xv[2] = angle_buffer[j];
        }
        //particles[i].xv[2] = pi_to_pi_base(xv2 + Vn*dt*tscheb_sin(Gn)/WHEELBASE);
    }
}
#endif

double predict_update_active_flops(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    int iwp_reset = *iwp;
    double G_reset = *G;
    double _xtrue[3] = {xtrue[0], xtrue[1], xtrue[2]};
    double flop_count = compute_steering_active_flops(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
    flop_count += predict_true_base_flops(V, *G, WHEELBASE, dt, xtrue);
    double VnGn[2];
    flop_count +=add_control_noise_base_flops(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn);
    //add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn);
    // flop_count += add_control_noise_base_flops(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    if (SWITCH_PREDICT_NOISE) {
        flop_count+=N*(2*tp.fastrand + 2* mmadd_2x2_flops(Q,Q,Q));
    }
    flop_count+= N* (5*tp.mul + 4* tp.add+
                    2*tp.sin /*tscheb_sin*/ + 1* tp.cos /*tscheb_cos*/+ 
                    pi_to_pi_active_flops(3.0) /*simple pi_to_pi*/);
    *iwp = iwp_reset;
    *G = G_reset;
    copy(_xtrue, 3, xtrue);
    return flop_count;
}

double predict_update_active_memory(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    int iwp_reset = *iwp;
    double G_reset = *G;
    double memory_moved = compute_steering_active_memory(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
    memory_moved += predict_true_base_memory(V, *G, WHEELBASE, dt, xtrue);
    // memory_moved += add_control_noise_base_memory(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    //add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn);
    memory_moved += 2*3*N; //We load xv and write to xv //The angle_buffer[4] is supposed to stay in cache, so I didnt count that
    *iwp = iwp_reset;
    *G = G_reset;
    return memory_moved;
}

#ifdef __AVX2__
void predict_update_fast(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
                        
    compute_steering_base(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
            if ( *iwp == -1 && NUMBER_LOOPS > 1 ) {
                *iwp = 0;
                NUMBER_LOOPS--;
            }
    predict_true_base(V, *G, WHEELBASE, dt, xtrue);

    // add process noise
    double VnGn[2] __attribute__ ((aligned(32)));

    add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    // Predict step
    double Vn, Gn, Vn1, Gn1, Vn2, Gn2, Vn3, Gn3;
    __m256d Vns, VndtWB, Vndt;
    __m256d Gns, Gn_theta, Gnsin;

    __m256d angles;

    __m256d dtv = _mm256_set1_pd(dt);
    __m256d WBv = _mm256_set1_pd(WHEELBASE);
    __m256d dtWBv = _mm256_div_pd(dtv, WBv);

    __m256i load_mask = _mm256_set_epi64x(9,6,3,0);

    __m256d thetas, xs, ys, xv1, xv2, xv3, xv, Gn_theta_sin, Gn_theta_cos;
    Vns = _mm256_set1_pd(VnGn[0]);
    Gns = _mm256_set1_pd(VnGn[1]);
    double angle_buffer[4] __attribute__ ((aligned(32)));

    double S[4]  __attribute__ ((aligned(32)));

    llt_2x2(Q,S);

    __m256d SMat = _mm256_load_pd(S);
    __m256d VG = _mm256_set_pd(VnGn[1], VnGn[0], VnGn[1], VnGn[0]);


    for (size_t i = 0; i < N; i+=4) {
        if (SWITCH_PREDICT_NOISE == 1) {
            //Vector2d noise;
            //multivariate_gauss_base(VnGn,Q,noise);	
            __m256d rand_vec1 = fill_rand_avx(-1.0,1.0);
            __m256d rand_vec2 = fill_rand_avx(-1.0,1.0);

            __m256d mul_gauss1 = _mmTadd_2x2_avx_v2(rand_vec1, SMat, VG);
            __m256d mul_gauss2 = _mmTadd_2x2_avx_v2(rand_vec2, SMat, VG);

            // S*rand + VnGn

            Vns = _mm256_shuffle_pd(mul_gauss1, mul_gauss2,0b0000);
            Gns = _mm256_shuffle_pd(mul_gauss1, mul_gauss2,0b1111);
        }
        xs = _mm256_i64gather_pd(particles[i].xv, load_mask, 8);
        ys = _mm256_i64gather_pd(particles[i].xv + 1, load_mask, 8);

        thetas = _mm256_i64gather_pd(particles[i].xv +2, load_mask, 8);

        Vndt = _mm256_mul_pd(Vns, dtv);
        VndtWB = _mm256_mul_pd(Vns, dtWBv);
        Gn_theta = _mm256_add_pd(Gns, thetas);
        Gn_theta_cos = tscheb_cos_avx(Gn_theta);
        Gn_theta_sin = tscheb_sin_avx(Gn_theta);
        Gnsin = tscheb_sin_avx(Gns);

#ifdef __FMA__
        xs = _mm256_fmadd_pd(Vndt, Gn_theta_cos, xs);
#else
        xs = _mm256_add_pd( _mm256_mul_pd(Vndt, Gn_theta_cos), xs);
#endif
#ifdef __FMA__
        ys = _mm256_fmadd_pd(Vndt, Gn_theta_sin, ys);
#else
        ys = _mm256_add_pd( _mm256_mul_pd(Vndt, Gn_theta_sin), ys);
#endif

        
#ifdef __FMA__
        angles = _mm256_fmadd_pd(VndtWB, Gnsin, thetas);
#else
        angles = _mm256_add_pd( _mm256_mul_pd(VndtWB, Gnsin), thetas);
#endif
        
        xv1 = _mm256_shuffle_pd(xs,ys,0b0000); // x1, y1, x3, y3
        xv2 = _mm256_shuffle_pd(xs,ys,0b1111); // x2, y2, x4, y4

        _mm256_store2_m128d(particles[i+2].xv, particles[i].xv, xv1);
        _mm256_store2_m128d(particles[i+3].xv, particles[i+1].xv, xv2);
        
        angles = simple_pi_to_pi_avx(angles);
        _mm256_store_pd(angle_buffer, angles);
        //double xv2 = particles[i].xv[2];
        //particles[i].xv[0] += Vn*dt*tscheb_cos(Gn + xv2);
        //particles[i].xv[1] += Vn*dt*tscheb_sin(Gn + xv2); 

        for (int j = 0; j<4; j++) {
            particles[i+j].xv[2] = angle_buffer[j];
        }
        //particles[i].xv[2] = pi_to_pi_base(xv2 + Vn*dt*tscheb_sin(Gn)/WHEELBASE);
    }
}
#endif


#ifdef __AVX2__
void predict_update_fast_normal_rand(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
                        
    compute_steering_base(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
            if ( *iwp == -1 && NUMBER_LOOPS > 1 ) {
                *iwp = 0;
                NUMBER_LOOPS--;
            }
    predict_true_base(V, *G, WHEELBASE, dt, xtrue);

    // add process noise
    double VnGn[2] __attribute__ ((aligned(32)));

    add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    // Predict step
    double Vn, Gn, Vn1, Gn1, Vn2, Gn2, Vn3, Gn3;
    __m256d Vns, VndtWB, Vndt;
    __m256d Gns, Gn_theta, Gnsin;

    __m256d angles;

    __m256d dtv = _mm256_set1_pd(dt);
    __m256d WBv = _mm256_set1_pd(WHEELBASE);
    __m256d dtWBv = _mm256_div_pd(dtv, WBv);

    __m256i load_mask = _mm256_set_epi64x(9,6,3,0);

    __m256d thetas, xs, ys, xv1, xv2, xv3, xv, Gn_theta_sin, Gn_theta_cos;
    Vns = _mm256_set1_pd(VnGn[0]);
    Gns = _mm256_set1_pd(VnGn[1]);
    double angle_buffer[4] __attribute__ ((aligned(32)));

    double S[4]  __attribute__ ((aligned(32)));

    llt_2x2(Q,S);
    __m256d SMat = _mm256_load_pd(S);
    __m256d VG = _mm256_set_pd(VnGn[1], VnGn[0], VnGn[1], VnGn[0]);

    double rng[8]  __attribute__ ((aligned(32)));
    for (size_t i = 0; i < N; i+=4) {
        if (SWITCH_PREDICT_NOISE == 1) {
            //BEGIN multivariate_gauss
            fill_rand(rng, 8, -1.0,1.0);


            __m256d rand_vec1 = _mm256_i64gather_pd(rng, _mm256_set_epi64x(5,4,1,0), 8);
            __m256d rand_vec2 = _mm256_i64gather_pd(rng,  _mm256_set_epi64x(7,6,3,2), 8);
            
            __m256d mul_gauss1 = _mmTadd_2x2_avx_v2(rand_vec1, SMat, VG);
            __m256d mul_gauss2 = _mmTadd_2x2_avx_v2(rand_vec2, SMat, VG);

            Vns = _mm256_shuffle_pd(mul_gauss1, mul_gauss2,0b0000);
            Gns = _mm256_shuffle_pd(mul_gauss1, mul_gauss2,0b1111);
            //END multivariate_gauss
        }

        xs = _mm256_i64gather_pd(particles[i].xv, load_mask, 8);
        ys = _mm256_i64gather_pd(particles[i].xv + 1, load_mask, 8);

        thetas = _mm256_i64gather_pd(particles[i].xv +2, load_mask, 8);

        Vndt = _mm256_mul_pd(Vns, dtv);
        VndtWB = _mm256_mul_pd(Vns, dtWBv);
        Gn_theta = _mm256_add_pd(Gns, thetas);
        Gn_theta_cos = tscheb_cos_avx(Gn_theta);
        Gn_theta_sin = tscheb_sin_avx(Gn_theta);
        Gnsin = tscheb_sin_avx(Gns);

#ifdef __FMA__
        xs = _mm256_fmadd_pd(Vndt, Gn_theta_cos, xs);
#else
        xs = _mm256_add_pd( _mm256_mul_pd(Vndt, Gn_theta_cos), xs);
#endif
#ifdef __FMA__
        ys = _mm256_fmadd_pd(Vndt, Gn_theta_sin, ys);
#else
        ys = _mm256_add_pd( _mm256_mul_pd(Vndt, Gn_theta_sin), ys);
#endif

        
#ifdef __FMA__
        angles = _mm256_fmadd_pd(VndtWB, Gnsin, thetas);
#else
        angles = _mm256_add_pd( _mm256_mul_pd(VndtWB, Gnsin), thetas);
#endif
        
        xv1 = _mm256_shuffle_pd(xs,ys,0b0000); // x1, y1, x3, y3
        xv2 = _mm256_shuffle_pd(xs,ys,0b1111); // x2, y2, x4, y4

        _mm256_store2_m128d(particles[i+2].xv, particles[i].xv, xv1);
        _mm256_store2_m128d(particles[i+3].xv, particles[i+1].xv, xv2);
        

        _mm256_store_pd(angle_buffer, angles);
        //double xv2 = particles[i].xv[2];
        //particles[i].xv[0] += Vn*dt*tscheb_cos(Gn + xv2);
        //particles[i].xv[1] += Vn*dt*tscheb_sin(Gn + xv2); 

        for (int j = 0; j<4; j++) {
            particles[i+j].xv[2] = pi_to_pi_base(angle_buffer[j]);
        }

        //particles[i].xv[2] = pi_to_pi_base(xv2 + Vn*dt*tscheb_sin(Gn)/WHEELBASE);
    }

}
#endif

#ifdef __AVX2__
void predict_update_fast_scalar_pipi(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
                        
    compute_steering_base(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
            if ( *iwp == -1 && NUMBER_LOOPS > 1 ) {
                *iwp = 0;
                NUMBER_LOOPS--;
            }
    predict_true_base(V, *G, WHEELBASE, dt, xtrue);

    // add process noise
    double VnGn[2] __attribute__ ((aligned(32)));

    add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    // Predict step
    double Vn, Gn, Vn1, Gn1, Vn2, Gn2, Vn3, Gn3;
    __m256d Vns, VndtWB, Vndt;
    __m256d Gns, Gn_theta, Gnsin;

    __m256d angles;

    __m256d dtv = _mm256_set1_pd(dt);
    __m256d WBv = _mm256_set1_pd(WHEELBASE);
    __m256d dtWBv = _mm256_div_pd(dtv, WBv);

    __m256i load_mask = _mm256_set_epi64x(9,6,3,0);

    __m256d thetas, xs, ys, xv1, xv2, xv3, xv, Gn_theta_sin, Gn_theta_cos;
    Vns = _mm256_set1_pd(VnGn[0]);
    Gns = _mm256_set1_pd(VnGn[1]);
    double angle_buffer[4] __attribute__ ((aligned(32)));

    double S[4]  __attribute__ ((aligned(32)));

    llt_2x2(Q,S);
    __m256d SMat = _mm256_load_pd(S);
    __m256d VG = _mm256_set_pd(VnGn[1], VnGn[0], VnGn[1], VnGn[0]);


    for (size_t i = 0; i < N; i+=4) {
        if (SWITCH_PREDICT_NOISE == 1) {
            //Vector2d noise;
            //multivariate_gauss_base(VnGn,Q,noise);	
            __m256d rand_vec1 = fill_rand_avx(-1.0,1.0);
            __m256d rand_vec2 = fill_rand_avx(-1.0,1.0);

            __m256d mul_gauss1 = _mmTadd_2x2_avx_v2(rand_vec1, SMat, VG);
            __m256d mul_gauss2 = _mmTadd_2x2_avx_v2(rand_vec2, SMat, VG);

            // S*rand + VnGn

            Vns = _mm256_shuffle_pd(mul_gauss1, mul_gauss2,0b0000);
            Gns = _mm256_shuffle_pd(mul_gauss1, mul_gauss2,0b1111);
        }
        xs = _mm256_i64gather_pd(particles[i].xv, load_mask, 8);
        ys = _mm256_i64gather_pd(particles[i].xv + 1, load_mask, 8);

        thetas = _mm256_i64gather_pd(particles[i].xv +2, load_mask, 8);

        Vndt = _mm256_mul_pd(Vns, dtv);
        VndtWB = _mm256_mul_pd(Vns, dtWBv);
        Gn_theta = _mm256_add_pd(Gns, thetas);
        Gn_theta_cos = tscheb_cos_avx(Gn_theta);
        Gn_theta_sin = tscheb_sin_avx(Gn_theta);
        Gnsin = tscheb_sin_avx(Gns);

#ifdef __FMA__
        xs = _mm256_fmadd_pd(Vndt, Gn_theta_cos, xs);
#else
        xs = _mm256_add_pd( _mm256_mul_pd(Vndt, Gn_theta_cos), xs);
#endif
#ifdef __FMA__
        ys = _mm256_fmadd_pd(Vndt, Gn_theta_sin, ys);
#else
        ys = _mm256_add_pd( _mm256_mul_pd(Vndt, Gn_theta_sin), ys);
#endif

        
#ifdef __FMA__
        angles = _mm256_fmadd_pd(VndtWB, Gnsin, thetas);
#else
        angles = _mm256_add_pd( _mm256_mul_pd(VndtWB, Gnsin), thetas);
#endif
        
        xv1 = _mm256_shuffle_pd(xs,ys,0b0000); // x1, y1, x3, y3
        xv2 = _mm256_shuffle_pd(xs,ys,0b1111); // x2, y2, x4, y4

        _mm256_store2_m128d(particles[i+2].xv, particles[i].xv, xv1);
        _mm256_store2_m128d(particles[i+3].xv, particles[i+1].xv, xv2);
        
        
        _mm256_store_pd(angle_buffer, angles);
        //double xv2 = particles[i].xv[2];
        //particles[i].xv[0] += Vn*dt*tscheb_cos(Gn + xv2);
        //particles[i].xv[1] += Vn*dt*tscheb_sin(Gn + xv2); 

        for (int j = 0; j<4; j++) {
            particles[i+j].xv[2] = pi_to_pi_base(angle_buffer[j]);
        }
        //particles[i].xv[2] = pi_to_pi_base(xv2 + Vn*dt*tscheb_sin(Gn)/WHEELBASE);
    }
}
#endif

#ifdef __AVX2__
void predict_update_simd(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    compute_steering_base(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
            if ( *iwp == -1 && NUMBER_LOOPS > 1 ) {
                *iwp = 0;
                NUMBER_LOOPS--;
            }
    predict_true_base(V, *G, WHEELBASE, dt, xtrue);

    // add process noise
    double VnGn[2] __attribute__ ((aligned(32)));

    add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    // Predict step
    double Vn, Gn, Vn1, Gn1, Vn2, Gn2, Vn3, Gn3;
    __m256d Vns, VndtWB, Vndt;
    __m256d Gns, Gn_theta, Gnsin;

    __m256d angles;

    __m256d dtv = _mm256_set1_pd(dt);
    __m256d WBv = _mm256_set1_pd(WHEELBASE);
    __m256d dtWBv = _mm256_div_pd(dtv, WBv);

    __m256i load_mask = _mm256_set_epi64x(9,6,3,0);

    __m256d thetas, xs, ys, xv1, xv2, xv3, xv, Gn_theta_sin, Gn_theta_cos;
    Vns = _mm256_set1_pd(VnGn[0]);
    Gns = _mm256_set1_pd(VnGn[1]);
    double angle_buffer[4] __attribute__ ((aligned(32)));
    for (size_t i = 0; i < N; i+=4) {
        if (SWITCH_PREDICT_NOISE == 1) {
            Vector2d noise;
            multivariate_gauss_base(VnGn,Q,noise);	
            Vn = noise[0];
            Gn = noise[1];
            Vector2d noise1;
            multivariate_gauss_base(VnGn,Q,noise1);	
            Vn1 = noise1[0];
            Gn1 = noise1[1];
            Vector2d noise2;
            multivariate_gauss_base(VnGn,Q,noise2);	
            Vn2 = noise2[0];
            Gn2 = noise2[1];
            Vector2d noise3;
            multivariate_gauss_base(VnGn,Q,noise3);	
            Vn3 = noise3[0];
            Gn3 = noise3[1];
            Vns = _mm256_set_pd(Vn3, Vn2, Vn1, Vn);
            Gns = _mm256_set_pd(Gn3, Gn2, Gn1, Gn);
        }
        xs = _mm256_i64gather_pd(particles[i].xv, load_mask, 8);
        ys = _mm256_i64gather_pd(particles[i].xv + 1, load_mask, 8);

        thetas = _mm256_i64gather_pd(particles[i].xv +2, load_mask, 8);

        Vndt = _mm256_mul_pd(Vns, dtv);
        VndtWB = _mm256_mul_pd(Vns, dtWBv);
        Gn_theta = _mm256_add_pd(Gns, thetas);
        Gn_theta_cos = tscheb_cos_avx(Gn_theta);
        Gn_theta_sin = tscheb_sin_avx(Gn_theta);
        Gnsin = tscheb_sin_avx(Gns);

#ifdef __FMA__
        xs = _mm256_fmadd_pd(Vndt, Gn_theta_cos, xs);
#else
        xs = _mm256_add_pd( _mm256_mul_pd(Vndt, Gn_theta_cos), xs);
#endif
#ifdef __FMA__
        ys = _mm256_fmadd_pd(Vndt, Gn_theta_sin, ys);
#else
        ys = _mm256_add_pd( _mm256_mul_pd(Vndt, Gn_theta_sin), ys);
#endif

        
#ifdef __FMA__
        angles = _mm256_fmadd_pd(VndtWB, Gnsin, thetas);
#else
        angles = _mm256_add_pd( _mm256_mul_pd(VndtWB, Gnsin), thetas);
#endif
        
        xv1 = _mm256_shuffle_pd(xs,ys,0b0000); // x1, y1, x3, y3
        xv2 = _mm256_shuffle_pd(xs,ys,0b1111); // x2, y2, x4, y4

        _mm256_store2_m128d(particles[i+2].xv, particles[i].xv, xv1);
        _mm256_store2_m128d(particles[i+3].xv, particles[i+1].xv, xv2);
        
        
        _mm256_store_pd(angle_buffer, angles);
        //double xv2 = particles[i].xv[2];
        //particles[i].xv[0] += Vn*dt*tscheb_cos(Gn + xv2);
        //particles[i].xv[1] += Vn*dt*tscheb_sin(Gn + xv2); 

        for (int j = 0; j<4; j++) {
            particles[i+j].xv[2] = pi_to_pi_base(angle_buffer[j]);
        }
        //particles[i].xv[2] = pi_to_pi_base(xv2 + Vn*dt*tscheb_sin(Gn)/WHEELBASE);
    }
}
#endif

void predict_update_sine(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    compute_steering_base(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
            if ( *iwp == -1 && NUMBER_LOOPS > 1 ) {
                *iwp = 0;
                NUMBER_LOOPS--;
            }
    predict_true_base(V, *G, WHEELBASE, dt, xtrue);

    // add process noise
    double VnGn[2];
    add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    // Predict step
    double Vn = V;
    double Gn = *G;
    for (size_t i = 0; i < N; i++) {
        if (SWITCH_PREDICT_NOISE == 1) {
            Vector2d noise;
            multivariate_gauss_base(VnGn,Q,noise);	
            Vn = noise[0];
            Gn = noise[1];
        }	
    
        double xv2 = particles[i].xv[2];
        particles[i].xv[0] += Vn*dt*tscheb_cos(Gn + xv2);
        particles[i].xv[1] += Vn*dt*tscheb_sin(Gn + xv2); 
        particles[i].xv[2] = pi_to_pi_base(xv2 + Vn*dt*tscheb_sin(Gn)/WHEELBASE);
    }
}

void predict_update_old(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    compute_steering(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
            if ( *iwp == -1 && NUMBER_LOOPS > 1 ) {
                *iwp = 0;
                NUMBER_LOOPS--;
            }
    predict_true(V, *G, WHEELBASE, dt, xtrue);

    // add process noise
    double VnGn[2];
    add_control_noise(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    // Predict step
    double Vn = V;
    double Gn = *G;

    for (size_t i = 0; i < N; i++) {
        if (SWITCH_PREDICT_NOISE == 1) {
            Vector2d noise;
            multivariate_gauss_base(VnGn,Q,noise);	
            Vn = noise[0];
            Gn = noise[1];
        }	
    
        double xv2 = particles[i].xv[2];
        particles[i].xv[0] += Vn*dt*cos(Gn + xv2);
        particles[i].xv[1] += Vn*dt*sin(Gn + xv2); 
        particles[i].xv[2] = pi_to_pi_base(xv2 + Vn*dt*sin(Gn)/WHEELBASE);
    }
}
