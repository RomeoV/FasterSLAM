#include "add_control_noise.h"
#include <iostream>
#include "linalg.h"
#include <math.h>

void add_control_noise(double V, double G, double* Q, int addnoise, double* VnGn)  {
	add_control_noise_active(V,G,Q,addnoise,VnGn);
}

void add_control_noise_base(double V, double G, double* Q, int addnoise, double* VnGn) 
{
	if (addnoise == 1) {
		Vector2d A; //seems we don't use malloc //malloc(2 * sizeof(double))
		A[0] = V;
		A[1] = G;
		Vector2d result; // need allocation? double *C = malloc(2 * sizeof(double));
		multivariate_gauss(A, Q, result);
		VnGn[0] = result[0];
		VnGn[1] = result[1];
	}
}

// Work / Memory instrumenting
void add_control_noise_active(double V, double G, double* Q, int addnoise, double* VnGn){
	if (addnoise == 1) {
		Vector2d result; // need allocation? double *C = malloc(2 * sizeof(double));
		
    // multivariate_gauss_active(A, Q, result); x P result
    /* inlining multivariate_gauss_active */
    double S[4]; //! 2x2 matrix, lower triangular cholesky factor
    llt_2x2(Q, S); //! P = S * S^T

    double X[2]; //! 2-vector
    fill_rand(X, 2, -1.0, 1.0);
    /* end inlining multivariate_gauss_active */

    //! result = S*X + x
    result[0] = V;
    result[1] = G;
    mvadd_2x2(S, X, result);

		VnGn[0] = result[0];
		VnGn[1] = result[1];
	}
}

double add_control_noise_base_flops(double V, double G, double* Q, int addnoise, double* VnGn){
  if(addnoise == 0){
	return 0;
  }

  Vector2d A;
  Vector2d result;
  double flop_count = multivariate_gauss_base_flops(A, Q, result);
  return flop_count;
}


double add_control_noise_base_memory(double V, double G, double* Q, int addnoise, double* VnGn){
  if(addnoise == 0){
	return 0;
  }

  Vector2d A;
  Vector2d result;
  double memory_called = multivariate_gauss_base_memory(A, Q, result);
  double memory_read_count = 2;
  double memory_written_count = 2 * 4;
  return memory_called + memory_read_count + memory_written_count;               
}

double add_control_noise_active_flops(double V, double G, double* Q, int addnoise, double* VnGn){
  if(addnoise == 0){
	return 0;
  }

  Vector2d A;
  Vector2d result;
  double flop_count = multivariate_gauss_active_flops(A, Q, result);
  return flop_count;
}

double add_control_noise_active_memory(double V, double G, double* Q, int addnoise, double* VnGn){
  if(addnoise == 0){
	return 0;
  }

  Vector2d A;
  Vector2d result;
  double memory_called = multivariate_gauss_active_memory(A, Q, result);
  double memory_read_count = 2;
  double memory_written_count = 2 * 4;
  return memory_called + memory_read_count + memory_written_count;  
}
