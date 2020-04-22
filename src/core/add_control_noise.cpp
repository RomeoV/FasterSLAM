#include "add_control_noise.h"

void add_control_noise(double V, double G, double* Q, int addnoise, double* VnGn) 
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
