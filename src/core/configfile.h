#pragma once

//******************
// Global Variables
//******************

extern double V;
extern double MAXG;
extern double RATEG;
extern double WHEELBASE;
extern double DT_CONTROLS;

extern double sigmaV;
extern double sigmaG;

extern double Q[2][2];

extern double MAX_RANGE;
extern double DT_OBSERVE;

extern double sigmaR;
extern double sigmaB;

extern double R[2][2];

extern double AT_WAYPOINT;
extern int NUMBER_LOOPS;

extern unsigned int NPARTICLES;
extern double NEFFECTIVE;

extern int SWITCH_CONTROL_NOISE;
extern int SWITCH_SENSOR_NOISE;
extern int SWITCH_INFLATE_NOISE;
extern int SWITCH_PREDICT_NOISE;

extern int SWITCH_SAMPLE_PROPOSAL;
extern int SWITCH_HEADING_KNOWN;
extern int SWITCH_RESAMPLE;
extern int SWITCH_PROFILE;
extern int SWITCH_SEED_RANDOM;