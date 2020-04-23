#include "predict_true.h"
#include "pi_to_pi.h"
#include <math.h>

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

void predict_true(const Vehicle* vehicle, const double steering_angle, const double dt, Vector3d* xv) {
    (*xv)[0] = vehicle->xtrue[0] + vehicle->V * dt * cos(vehicle->xtrue[2] + steering_angle);
    (*xv)[1] = vehicle->xtrue[1] + vehicle->V * dt * sin(vehicle->xtrue[2] + steering_angle);
    (*xv)[2] = pi_to_pi(vehicle->xtrue[2] + vehicle->V*dt*sin(steering_angle)/(-vehicle->veh[0][1]));		
}
