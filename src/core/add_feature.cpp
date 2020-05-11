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
    add_feature_base(particle,z, N_z, R);
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
  // double* xf = (double*)malloc(2 * N_z * sizeof(double));
  // double* Pf = (double*)malloc(4 * N_z * sizeof(double));
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

  // free(xf);
  // free(Pf);
}
