#ifndef ACOUSTICS_CONST_RP_H
#define ACOUSTICS_CONST_RP_H

#include <math.h>
#include <stdio.h>

typedef double real;

struct rp_params
{
    // Bulk modules (k)
    real bulk;
    // Density
    real rho;
    // Impedence
    real Z;
    // Speed of sound
    real c
}

// These values are all constant for this Riemann solver so calculate them
// here instead of in the Riemann solver
static const rp_params acoustics_rp_params = 
{ 
    bulk = 1.0;
    rho = 1.0;
    c = sqrt(bulk / rho);
    Z = c * rho ;
}

// Point-wise constant acoustics Riemann solver
inline
void acoustics_const_rp(const real* q_left, const real* q_right,
                        const real* aux_left, const real* aux_right,
                        rp_params* params,
                        real* amdq, real* apdq, real *wave, real *s)
{
    // Local physical constants
    real[2] alpha;
    
    // Wave strengths
    real alpha[0] = ( q_left[0] - q_right[0] + 
                        params.z * (q_right[1] - q_left[1])) / (2.0 * params.z);
    real alpha[1] = ( q_right[0] - q_left[0] +
                        params.z * (q_right[1] - q_left[1])) / (2.0 * params.z);
    
    // Wave speeds
    s[0] = -params.c, s[1] = params.c;
 
    // Set wave vectors
    wave[0,0] = -alpha[0] * params.z; 
    wave[0,1] =  alpha[0];
    wave[0,2] =  0.0;
    
    wave[1,0] = alpha[1] * params.z;
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