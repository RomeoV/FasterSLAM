#include "typedefs.h"

#include "particle.h"

/*! Thrun03 Eq. 61
    Compute particle weight for sampling.
    Uses compute_jacobians.
    @param[out] particle Particle whose weight is calculated
    @param[in]  z        vector of map features, calculated by data_associate_known
    @param[in]  idf      vector of map indices, calculated by data_associate_known, used for jacobians 
    @param[in]  R        matrix of observation noises, metres and radians
 */
double compute_weight(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R);