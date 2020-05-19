#include "add_control_noise.h"


void add_control_noise(double V, double G, double* Q, int addnoise, double* VnGn)  {
	add_control_noise_base(V,G,Q,addnoise,VnGn);
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
void add_control_noise_active(double V, double G, double* Q, int addnoise, double* VnGn);

double add_control_noise_base_flops(double V, double G, double* Q, int addnoise, double* VnGn);

double add_control_noise_base_memory(double V, double G, double* Q, int addnoise, double* VnGn);

double add_control_noise_active_flops(double V, double G, double* Q, int addnoise, double* VnGn);

double add_control_noise_active_memory(double V, double G, double* Q, int addnoise, double* VnGn);
