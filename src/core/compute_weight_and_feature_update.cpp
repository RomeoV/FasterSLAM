
#include "compute_weight_and_feature_update.h"
#include "compute_weight.h"
#include "feature_update.h"

#include "linalg.h"
#include <math.h>
#include "compute_jacobians.h"


void compute_weight_and_feature_update(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R) {
    compute_weight_and_feature_update_base(particle, z, idf, N_idf, R);
}

void compute_weight_and_feature_update_base(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R) {
    Vector2d zp[N_idf];
    Matrix23d Hv[N_idf];
    Matrix2d Hf[N_idf];
    Matrix2d Sf[N_idf];
    *(particle->w) *= compute_weight_base(particle, z,N_idf, idf, R, zp, Hv, Hf, Sf);
    feature_update_base(particle, z, idf,N_idf, R, zp, Hv, Hf, Sf);
}

void compute_weight_and_feature_update_active(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R){
    Vector2d zp[N_idf] __attribute__((aligned(32)));
    Matrix23d Hv[N_idf] __attribute__((aligned(32)));
    Matrix2d Hf[N_idf] __attribute__((aligned(32)));
    Matrix2d Sf[N_idf] __attribute__((aligned(32)));

    // process each feature, incrementally refine proposal distribution
    compute_jacobians_fast(particle, idf, N_idf, R, zp, Hv, Hf, Sf);

    Vector2d v[N_idf];
    for (size_t j = 0; j < N_idf; j++) {
        Vector2d v_j;
        sub(z[j], zp[j], 2, v_j);  // v_j = z[j] - zp[j]
        v_j[1] = pi_to_pi(v_j[1]);
        copy(v_j, 2, v[j]);  // v[j] = v_j
    }

    double w = 1.0;

    double den, num;
    // this can probably be done alot faster without this loop.....
    
    for (size_t i = 0; i < N_idf; i++) {
        Matrix2d S, ST, S_inv;
        Vector2d S_inv_v;
        double vT_S_inv_v;
        // Eq. 61 in Thrun03g
        
        copy(Sf[i], 4, S);
        transpose(S, 2, 2, ST);
        inv_2x2(S, S_inv);

        
        mv_2x2(S_inv, v[i], S_inv_v);
        mul(v[i], S_inv_v, 1, 2, 1, &vT_S_inv_v); // TODO in linalg

        

        den = 2 * M_PI * sqrt(determinant_2x2(S));
        num = exp(-0.5 * vT_S_inv_v);
        w *= (double)num / (double)den;
    }

    *(particle->w)*=w;

    ///////////////////FEATURE_UPDATE//////////////////////////////////////////

      // Having selected a new pose from the proposal distribution, this pose is
    // assumed perfect and each feature update maybe computed independently and
    // without pose uncertainty

    // double* xf = (double*)malloc(2 * N_idf * sizeof(double));
    // double* Pf = (double*)malloc(4 * N_idf * sizeof(double));
    Vector2d xf[N_idf];
    Matrix2d Pf[N_idf];

    for (size_t i = 0; i < N_idf; i++) {
        copy(particle->xf + (2 * idf[i]), 2, xf[i]);  // means
        copy(particle->Pf + (4 * idf[i]), 4, Pf[i]);  // covariances
    }

    Vector2d feat_diff[N_idf];  // difference btw feature prediciton and
                                // measurement (used to update mean)
    for (int i = 0; i < N_idf; i++) {
        sub(z[i], zp[i], 2, feat_diff[i]);
        feat_diff[i][1] = pi_to_pi_base(feat_diff[i][1]);
    }

    for (int i = 0; i < N_idf; i++) {
        KF_cholesky_update_base(xf[i], Pf[i], 
                        feat_diff[i], R, 
                        Hf[i]);
    }

    for (size_t i = 0; i < N_idf; i++) {
        set_xfi(particle, xf[i], idf[i]);
        set_Pfi(particle, Pf[i], idf[i]);
    }
}