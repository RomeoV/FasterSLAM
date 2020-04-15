#include "configfile.h"

#define pi 3.14159265

// Configuration file
//Permits various adjustments to parameters of the FastSLAM algorithm.
// See fastslam_sim.h for more information

// control parameters
double V= 3.0; // m/s
double MAXG= 30*pi / 180; // radians, maximum steering angle (-MAXG < g < MAXG)
double RATEG= 20*pi / 180; // rad/s, maximum rate of change in steer angle
double WHEELBASE= 4.; // metres, vehicle wheel-base
double DT_CONTROLS= 0.025; // seconds, time interval between control signals

// control noises
double sigmaV= 0.3; // m/s
double sigmaG= (3.0 * pi / 180); // radians

double Q[2][2]= {
    {0, 0},
    {0, 0}
};

// observation parameters
double MAX_RANGE= 30.0; // metres
double DT_OBSERVE= 8 * 0.025; //DT_CONTROLS; // seconds, time interval between observations

// observation noises
double sigmaR= 0.1; // metres
double sigmaB= (1.0 * pi / 180); // radians

double R[2][2] = {
    {0, 0},
    {0, 0}
};

// waypoint proximity
double AT_WAYPOINT= 1.0; // metres, distance from current waypoint at which to switch to next waypoint
int NUMBER_LOOPS= 2; // number of loops through the waypoint list

// resampling
unsigned int NPARTICLES= 100; 
double NEFFECTIVE= 0.75* 100; //NPARTICLES; // minimum number of effective particles before resampling

// switches
int SWITCH_CONTROL_NOISE= 1;
int SWITCH_SENSOR_NOISE = 1;
int SWITCH_INFLATE_NOISE= 0;
int SWITCH_PREDICT_NOISE = 1;   // sample noise from predict (usually 1 for fastslam1.0 and 0 for fastslam2.0)
int SWITCH_SAMPLE_PROPOSAL = 1; // sample from proposal (no effect on fastslam1.0 and usually 1 for fastslam2.0)
int SWITCH_HEADING_KNOWN= 0;
int SWITCH_RESAMPLE= 1; 
int SWITCH_PROFILE= 1;
int SWITCH_SEED_RANDOM= 0; // if not 0, seed the randn() with its value at beginning of simulation (for repeatability)

