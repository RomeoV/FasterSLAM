#include "add_control_noise.h"


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
		multivariate_gauss_base(A, Q, result);
		VnGn[0] = result[0];
		VnGn[1] = result[1];
	}
}

// Work / Memory instrumenting
void add_control_noise_active(double V, double G, double* Q, int addnoise, double* VnGn){
	if (addnoise == 1) {
		Vector2d A; //seems we don't use malloc //malloc(2 * sizeof(double))
		A[0] = V;
		A[1] = G;
		Vector2d result; // need allocation? double *C = malloc(2 * sizeof(double));
		multivariate_gauss_active(A, Q, result);
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
