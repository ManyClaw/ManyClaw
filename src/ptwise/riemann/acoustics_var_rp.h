#ifndef ACOUSTICS_VAR_RP_H
#define ACOUSTICS_VAR_RP_H

#include <math.h>
#include <stdio.h>

typedef double real;

typedef struct rp_params
{} rp_params;

// Point-wise constant acoustics Riemann solver
// aux_* = [rho,bulk]
inline
void acoustics_var_rp(const real* q_left, const real* q_right,
                      const real* aux_left, const real* aux_right,
                      rp_params* params,
                      real* amdq, real* apdq, real *wave, real *s)
{
    // Local physical constants
    real c_left = sqrt(aux_left[1] / aux_left[0]);
    real c_right = sqrt(aux_right[1] / aux_right[0]);
    real Z_left = c_left * aux_left[0];
    real Z_right = c_right * aux_right[0];
    
    // Wave strengths
    real alpha[2];
    alpha[0] = ( q_left[0] - q_right[0] + Z_right * (q_right[1] - q_left[1]) ) 
                        / (Z_left + Z_right);
    alpha[1] = ( q_right[0] - q_left[0] + Z_left * (q_right[1] - q_left[1]) ) 
                        / (Z_left + Z_right);
    
    // Wave speeds
    s[0] = -c_left;
    s[1] =  c_right;
 
    // Set wave vectors
    wave[0,0] = -alpha[0] * Z_left; 
    wave[0,1] = alpha[0];
    wave[0,2] = 0.0;
    
    wave[1,0] = alpha[1] * Z_right;
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

#endif // ACOUSTICS_VAR_RP_H