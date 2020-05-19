#include "compute_weight.h"

#include <math.h>

#include "compute_jacobians.h"
#include "linalg.h"
#include "pi_to_pi.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation
 * ToDo: Unit tests
 ****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: TBD
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/

double compute_weight(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]) {
    return compute_weight_base(particle, z, N_z, idf, R, zp, Hv, Hf, Sf);
}

double compute_weight_base(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]) {

  // process each feature, incrementally refine proposal distribution
  compute_jacobians_base(particle, idf, N_z, R, zp, Hv, Hf, Sf);

  Vector2d v[N_z];
  for (size_t j = 0; j < N_z; j++) {
    Vector2d v_j;
    sub(z[j], zp[j], 2, v_j);  // v_j = z[j] - zp[j]
    v_j[1] = pi_to_pi_base(v_j[1]);
    copy(v_j, 2, v[j]);  // v[j] = v_j
  }

  double w = 1.0;

  double den, num;
  // this can probably be done alot faster without this loop.....
  
  for (size_t i = 0; i < N_z; i++) {
    Matrix2d S, ST, S_inv;
    Vector2d S_inv_v;
    double vT_S_inv_v;
    // Eq. 61 in Thrun03g
    
    copy(Sf[i], 4, S);
    transpose(S, 2, 2, ST);
    inv_2x2(S, S_inv);
    mul(S_inv, v[i], 2, 2, 1, S_inv_v);
    mul(v[i], S_inv_v, 1, 2, 1, &vT_S_inv_v);

    den = 2 * M_PI * sqrt(determinant_2x2(S));
    num = exp(-0.5 * vT_S_inv_v);
    w *= (double)num / (double)den;
  }

  return w;
}

double compute_weight_active(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]) {

  // process each feature, incrementally refine proposal distribution

  Vector2d v[N_z];
  for (size_t j = 0; j < N_z; j++) {
    Vector2d v_j;
    sub(z[j], zp[j], 2, v_j);  // v_j = z[j] - zp[j]
    v_j[1] = pi_to_pi_active(v_j[1]);
    copy(v_j, 2, v[j]);  // v[j] = v_j
  }

  double w = 1.0;

  double den, num;
  // this can probably be done alot faster without this loop.....
  
  for (size_t i = 0; i < N_z; i++) {
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

  return w;
}

// Work / Memory instrumenting
double compute_weight_base_flops(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]){

  Vector2d v_j;
  Matrix2d S, ST, S_inv;
  Vector2d S_inv_v;
  copy(Sf[0], 4, S);
  transpose(S, 2, 2, ST);
  double vT_S_inv_v;

  double flop_count = compute_jacobians_base(particle, idf, N_z, R, zp, Hv, Hf, Sf) + N_z * (
      4 * tp.mul +
      1 * tp.div +
      1 * tp.sqrt + 
      1 * tp.exp + 
      sub_flops(z[0], zp[0], 2, v_j) + 
      pi_to_pi_base_flops(v_j[1]) + // TODO: depends on input
      inv_2x2_flops(S, S_inv) + 
      mul_flops(S_inv, v_j, 2, 2, 1, S_inv_v) +
      mul_flops(v_j, S_inv_v, 1, 2, 1, &vT_S_inv_v) + 
      determinant_2x2_flops(S)
    );
  return flop_count;
}


double compute_weight_base_memory(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]){

  Vector2d v[N_z];
  Vector2d v_j;
  Matrix2d S, ST, S_inv;
  Vector2d S_inv_v;
  double vT_S_inv_v;

  double memory_called = compute_jacobians_base(particle, idf, N_z, R, zp, Hv, Hf, Sf) +
  N_z * (
    sub_memory(z[0], zp[0], 2, v_j) + 
    pi_to_pi_base_memory(v_j[1]) +
    copy_memory(v_j, 2, v[j]) + // 2 * 2 + //copy_memory(v_j, 2, v[j]) + 
    copy_memory(Sf[i], 4, S) + // 2 * 4 + //copy_memory(Sf[i], 4, S) +
    transpose_memory(S, 2, 2, ST) +
    inv_2x2_memory(S, S_inv) +
    mul_memory(S_inv, v[0], 2, 2, 1, S_inv_v) +
    mul_memory(v[0], S_inv_v, 1, 2, 1, &vT_S_inv_v) +
    determinant_2x2_memory(S)
  );
  double memory_read_count = N_z * 7;
  double memory_written_count = N_z * 2;
  return memory_called + memory_read_count + memory_written_count;               
}

double compute_weight_active_flops(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]){

  Vector2d v_j;
  Matrix2d S, ST, S_inv;
  Vector2d S_inv_v;
  copy(Sf[0], 4, S);
  transpose(S, 2, 2, ST);
  double vT_S_inv_v;

  double flop_count = N_z * (
      4 * tp.mul +
      1 * tp.div +
      1 * tp.sqrt + 
      1 * tp.exp + 
      sub_flops(z[0], zp[0], 2, v_j) + 
      pi_to_pi_active_flops(v_j[1]) + // TODO: depends on input
      inv_2x2_flops(S, S_inv) + 
      mv_2x2_flops(S_inv, v_j, S_inv_v) + // TODO: where is this function
      mul_flops(v_j, S_inv_v, 1, 2, 1, &vT_S_inv_v) + 
      determinant_2x2_flops(S)
    );
  return flop_count;
}

double compute_weight_active_memory(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]){
  
  Vector2d v[N_z];
  Vector2d v_j;
  Matrix2d S, ST, S_inv;
  Vector2d S_inv_v;
  double vT_S_inv_v;

  double memory_called = N_z * (
    sub_memory(z[0], zp[0], 2, v_j) + 
    pi_to_pi_active_memory(v_j[1]) +
    copy_memory(v_j, 2, v[j]) + // 2 * 2 + //copy_memory(v_j, 2, v[j]) + 
    copy_memory(Sf[i], 4, S) + // 2 * 4 + //copy_memory(Sf[i], 4, S) +
    transpose_memory(S, 2, 2, ST) +
    inv_2x2_memory(S, S_inv) +
    mv_2x2_memory(S_inv, v[0], S_inv_v) + // TODO: where is this function
    mul_memory(v[0], S_inv_v, 1, 2, 1, &vT_S_inv_v) +
    determinant_2x2_memory(S)
  );
  double memory_read_count = N_z * 7;
  double memory_written_count = N_z * 2;
  return memory_called + memory_read_count + memory_written_count; 
}