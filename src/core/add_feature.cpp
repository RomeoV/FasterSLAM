#include "add_feature.h"

#include <assert.h>
#include <math.h>

#include "linalg.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
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
void add_feature(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R) {
    add_feature_base(particle,z, N_z, R);
}

void add_feature_base(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R) {
  // double* xf = (double*)malloc(2 * N_z * sizeof(double));
  // double* Pf = (double*)malloc(4 * N_z * sizeof(double));
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

  // free(xf);
  // free(Pf);
}
