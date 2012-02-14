#ifndef ACOUSTICS_CONST_RP_H
#define ACOUSTICS_CONST_RP_H

#include <math.h>
#include <stdio.h>

typedef double real;
const real BULK = 1.0;
const real RHO = 1.0;

typedef struct rp_params
{
    // Bulk modules (k)
    real bulk;
    // Density
    real rho;
    // Impedence
    real c;  // = sqrt(bulk / rho);
    // Speed of sound
    real Z; // = c * rho;
} rp_params;

// Point-wise constant acoustics Riemann solver
inline
void acoustics_const_rp(const real* q_left, const real* q_right,
                        const real* aux_left, const real* aux_right,
                        const rp_params* params,
                        real* amdq, real* apdq, real *wave, real *s)
{
    // Local physical constants
    real alpha[2];
    
    // Wave strengths
    alpha[0] = ( q_left[0] - q_right[0] + 
                        params->Z * (q_right[1] - q_left[1])) / (2.0 * params->Z);
    alpha[1] = ( q_right[0] - q_left[0] +
                        params->Z * (q_right[1] - q_left[1])) / (2.0 * params->Z);
    
    // Wave speeds
    s[0] = -params->c, s[1] = params->c;
 
    // Set wave vectors
    wave[0,0] = -alpha[0] * params->Z; 
    wave[0,1] =  alpha[0];
    wave[0,2] =  0.0;
    
    wave[1,0] = alpha[1] * params->Z;
    wave[1,1] = alpha[1];
    wave[1,2] = 0.0;
    
    // Grid edge fluctuations
    amdq[0] = s[0] * wave[0,0];
    amdq[1] = s[0] * wave[0,1];
    amdq[2] = s[0] * wave[0,2];  // This could just be set to zero
    
    apdq[0] = s[1] * wave[1,0];
    apdq[1] = s[1] * wave[1,1];
    apdq[2] = s[1] * wave[1,2];  // This could just be set to zero
}

#endif // ACOUSTICS_CONST_RP_H
