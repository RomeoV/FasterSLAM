#pragma once

typedef struct {
    double xtrue[3];
    double veh[2][3];
    double V;
    // double alpha; // initial steer angle, used to be called G in yglee
} Vehicle;
