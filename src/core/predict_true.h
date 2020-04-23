#pragma once
#include "typedefs.h"
#include "vehicle.h"

/*!
    Clips all angles in angle to range [-pi,pi]. [Simple to optimize]
    @param[in]  V   Velocity in m/s.
    @param[in]  G   Steering angle in radiants.
    @param[in]  WB  Wheelbase.
    @param[in]  dt  timestep.
    @param[out] xv  State vector (x,y,angle).
 */

void predict_true(const Vehicle* vehicle, const double steering_angle, const double dt, Vector3d* xv);
    

