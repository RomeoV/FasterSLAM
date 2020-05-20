#include "add_feature.h"

#include <assert.h>
#include <math.h>

#include "linalg.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/

void add_feature(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R) {
    add_feature_active(particle,z, N_z, R);
}

/*****************************************************************************
 * PERFORMANCE STATUS (N_z = nb of not before seen features)
 * Work: N_z*(4 adds + 2 muls + 2 sin/cos + 2 matmul 2x2x2 = 2*(8 muls + 4 adds))      + 2*N_z + 1 integer adds
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void add_feature_base(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R) {

  Vector2d xf[N_z];
  Matrix2d Pf[N_z];

  Vector3d xv;
  copy(particle->xv, 3, xv);

  double r, b, s, c;

  for (size_t i = 0; i < N_z; i++) {
    r = z[i][0];
    b = z[i][1];
    s = sin(xv[2] + b);
    c = cos(xv[2] + b);

    Vector2d measurement;
    measurement[0] = xv[0] + r * c;
    measurement[1] = xv[1] + r * s;
    copy(measurement, 2, xf[i]);  // xf[i,:] = measurement[0:2]

    Matrix2d Gz = {c, -r * s, s, r * c};
    Matrix2d MatResult_1;
    mul(Gz, R, 2, 2, 2, MatResult_1);
    Matrix2d Gz_T;
    transpose(Gz, 2, 2, Gz_T);
    Matrix2d MatResult_2;
    mul(MatResult_1, Gz_T, 2, 2, 2, MatResult_2);

    copy(MatResult_2, 2 * 2, Pf[i]);  // Pf[i,0:4] = Gz*R*Gz.T
  }

  size_t N_x = particle->Nfa;
  assert(particle->Nfa + N_z < particle->Nf);
  particle->Nfa += N_z;

  for (size_t i = 0; i < N_z; i++) {
    set_xfi(particle, xf[i], i + N_x);
    set_Pfi(particle, Pf[i], i + N_x);
  }

}

/*****************************************************************************
 * PERFORMANCE STATUS (N_z = nb of not before seen features)
 * Work: N_z*(4 adds + 2 muls + 2 sin/cos + 2 matmul 2x2x2 = 2*(8 muls + 4 adds))      + 2*N_z + 1 integer adds
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void add_feature_active(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R) {

  Vector2d xf[N_z] __attribute__((aligned(32)));
  Matrix2d Pf[N_z] __attribute__((aligned(32)));

  Vector3d xv __attribute__((aligned(32)));
  copy(particle->xv, 3, xv);

  double r, b, s, c;

  for (size_t i = 0; i < N_z; i++) {
    r = z[i][0];
    b = z[i][1];
    s = sin(xv[2] + b);
    c = cos(xv[2] + b);

    Vector2d measurement;
    measurement[0] = xv[0] + r * c;
    measurement[1] = xv[1] + r * s;
    copy(measurement, 2, xf[i]);  // xf[i,:] = measurement[0:2]

    Matrix2d Gz __attribute__((aligned(32))) = {c, -r * s, s, r * c} ;
    Matrix2d MatResult_1 __attribute__((aligned(32)));
    
#ifdef __AVX2__
    mm_2x2_avx_v1(Gz, R, MatResult_1);
    mmT_2x2_avx_v1(MatResult_1, Gz, Pf[i]);
#else
    mm_2x2(Gz, R, MatResult_1);
    mmT_2x2(MatResult_1, Gz, Pf[i]);
#endif

  }

  size_t N_x = particle->Nfa;
  assert(particle->Nfa + N_z < particle->Nf);
  particle->Nfa += N_z;

  for (size_t i = 0; i < N_z; i++) {
    set_xfi(particle, xf[i], i + N_x);
    set_Pfi(particle, Pf[i], i + N_x);
  }

}

// Work / Memory instrumenting
double add_feature_base_flops(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R){
  
  Matrix2d MatResult_1, MatResult_2, Gz, Gz_T;
  double flop_count = N_z * (
      tp.sin + tp.cos + 
      4*tp.add +  
      4*tp.mul +
      mul_flops(Gz, R, 2, 2, 2, MatResult_1) +
      mul_flops(MatResult_1, Gz_T, 2, 2, 2, MatResult_2)
      );
  
  return flop_count;
  }

double add_feature_base_memory(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R){

  Matrix2d MatResult_1, MatResult_2, Gz, Gz_T;
  Vector3d xv;
  Vector2d xf[N_z];
  Matrix2d Pf[N_z];
  Vector2d measurement;

  double memory_called = copy_memory(particle->xv, 3, xv) + N_z * (
    copy_memory(measurement, 2, xf[0]) +
    mul_memory(Gz, R, 2, 2, 2, MatResult_1) +
    transpose_memory(Gz, 2, 2, Gz_T) +
    mul_memory(MatResult_1, Gz_T, 2, 2, 2, MatResult_2) +
    copy_memory(MatResult_2, 2 * 2, Pf[0]) +
    2 * (2 + 1) + // set_xfi(particle, xf[i], i + N_x) + 
    2 * (2 + 1) // set_Pfi(particle, Pf[i], i + N_x)
  );
  double memory_read_count = N_z * 10;
  double memory_written_count = N_z * (2 * 2);
  return memory_called + memory_read_count + memory_written_count;
}

double add_feature_active_flops(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R){

  Matrix2d MatResult_1, Gz;
  Vector3d xv;
  Matrix2d Pf[N_z];

  double flop_count = N_z * (
      tp.sin + tp.cos + 
      4*tp.add +  
      4*tp.mul +
      /* different */
      // assuming we have AVX2
      mm_2x2_flops(Gz, R, MatResult_1) + // the same as mm_2x2_avx_v1_flops(Gz, R, MatResult_1) + 
      mm_2x2_flops(MatResult_1, Gz, Pf[0]) // the same as mmT_2x2_avx_v1_flops(MatResult_1, Gz, Pf[i])
      /* different end */
      );
  
  return flop_count;
  }

double add_feature_active_memory(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R){

  Matrix2d MatResult_1, Gz;
  Vector3d xv;
  Vector2d xf[N_z];
  Matrix2d Pf[N_z];
  Vector2d measurement;

  double memory_called = copy_memory(particle->xv, 3, xv) + N_z * (
    copy_memory(measurement, 2, xf[0]) +
    mm_2x2_flops(Gz, R, MatResult_1) + // the same as mm_2x2_avx_v1_memory(Gz, R, MatResult_1) + 
    mm_2x2_flops(MatResult_1, Gz, Pf[0]) + // the same as mmT_2x2_avx_v1_memory(MatResult_1, Gz, Pf[i]) +
    2 * (2 + 1) + // set_xfi(particle, xf[i], i + N_x) + 
    2 * (2 + 1) // set_Pfi(particle, Pf[i], i + N_x)
  );
  double memory_read_count = N_z * 10;
  double memory_written_count = N_z * (2 * 2);
  return memory_called + memory_read_count + memory_written_count;
  }
