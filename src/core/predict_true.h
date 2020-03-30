#pragma once
#include "typedefs.h"

/*!
    Clips all angles in angle to range [-pi,pi]. [Simple to optimize]
    @param[in]  V   Velocity in m/s.
    @param[in]  G   Steering angle in radiants.
    @param[in]  WB  Wheelbase.
    @param[in]  dt  timestep.
    @param[out] xv  State vector (x,y,angle).
 */

void predict_true(double V,double G,double WB,double dt, Vector3d xv);
    

