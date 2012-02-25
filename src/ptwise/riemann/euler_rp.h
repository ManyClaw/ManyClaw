#ifndef EULER_RP_H
#define EULER_RP_H

#include "../../many_claw.h"

struct euler_rp_aux_global_t
{
    real gamma;
};

static const euler_rp_aux_global_t euler_rp_aux_global = {1.0};
static const rp_grid_params euler_rp_grid_params = {2, 3, 0, 3};

// Evaluate the ideal gas law
inline real pressure(real rho, real momentum, real E)
{
    return 0.0;
}

// Implementation of the Roe solver for the Euler equations
inline
void euler_rp(const real* q_left, const real* q_right,
              const real* aux_left, const real* aux_right,
              const euler_rp_aux_global_t* aux_global,
              real* amdq, real* apdq, real* wave, real* s)
{
  const int num_states = euler_rp_grid_params.num_states;

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
  wave[0*num_states + 0] = 0.0;
  wave[0*num_states + 1] = 0.0;
  wave[0*num_states + 2] = 0.0;

  wave[1*num_states + 0] = 0.0;
  wave[1*num_states + 1] = 0.0;
  wave[1*num_states + 2] = 0.0;

  wave[2*num_states + 0] = 0.0;
  wave[2*num_states + 1] = 0.0;
  wave[2*num_states + 2] = 0.0;

  // Calculate fluctuations
  for (int eqn_index = 0; eqn_index < 3; eqn_index++)
  {
    amdq[eqn_index] = 0.0;
    apdq[eqn_index] = 0.0;
  }
  // Depending on the wave speed, amdq or apdq is modified
  for (int wave_index = 0; wave_index < 3; wave_index++)
  {
    if (s[wave_index] <= 0.0)
    {
      for (int eqn_index = 0; eqn_index < 3; eqn_index++)
      {
        amdq[eqn_index] += s[wave_index] * wave[wave_index + eqn_index*num_states];
      }
    }
    else
    {
      for (int eqn_index = 0; eqn_index < 3; eqn_index++)
      {
        apdq[eqn_index] += s[wave_index] * wave[wave_index + eqn_index*num_states];
      }
    }
  }
}

#endif
