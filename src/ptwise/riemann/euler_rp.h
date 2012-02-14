#ifndef EULER_RP_H
#define EULER_RP_H

#include <math.h>

typedef double real;

typedef struct rp_params
{
    real gamma;
} rp_params;

static const rp_params advection_rp_params = {1.4};

// Evaluate the ideal gas law
inline real pressure(real rho, real momentum, real E)
{
    return 0.0;
}

// Implementation of the Roe solver for the Euler equations
inline
void euler_rp(const real* q_left, const real* q_right,
              const real* aux_left, const real* aux_right,
              const rp_params* params,
              real* amdq, real* apdq, real* wave, real* s)
{   
    // Loop indices
    int eqn_index, wave_index;
    
    // Roe averages
    real rho_hat = 0.0;
    real rhou_hat = 0.0;
    real E_hat = 0.0;
    
    // Wave strengths
    real alpha[3];
    
    // Wave speeds
    s[0] = 0.0;
    s[1] = 0.0;
    s[2] = 0.0;
    
    // Waves
    wave[0,0] = 0.0;
    wave[0,1] = 0.0;
    wave[0,2] = 0.0;

    wave[1,0] = 0.0;
    wave[1,1] = 0.0;
    wave[1,2] = 0.0;
    
    wave[2,0] = 0.0;
    wave[2,1] = 0.0;
    wave[2,2] = 0.0;
    
    // Calculate fluctuations
    for (eqn_index = 0; eqn_index < 3; eqn_index++)
    {
        amdq[eqn_index] = 0.0;
        apdq[eqn_index] = 0.0;
    }
    // Depending on the wave speed, amdq or apdq is modified
    for (wave_index = 0; wave_index < 3; wave_index++)
    {
        if (s[wave_index] <= 0.0) 
        {
            for (eqn_index = 0; eqn_index < 3; eqn_index++)
            {
                amdq[eqn_index] += s[wave_index] * wave[wave_index,eqn_index];                
            }   
        }
        else 
        {
            for (eqn_index = 0; eqn_index < 3; eqn_index++)
            {
                apdq[eqn_index] += s[wave_index] * wave[wave_index,eqn_index];                
            }   
        }
    }
}

#endif
