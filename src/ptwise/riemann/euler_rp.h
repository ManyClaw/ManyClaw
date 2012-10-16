#ifndef EULER_RP_H
#define EULER_RP_H

#include "../../many_claw.h"

struct euler_rp_aux_global_t
{
    real gamma;
};

static const euler_rp_aux_global_t euler_rp_aux_global = {1.0};
static const rp_grid_params euler_rp_grid_params = {2, 4, 0, 4};

// Evaluate the ideal gas law
inline real pressure(real rho, real momentum, real E)
{
    return 0.0;
}

// Implementation of the Roe solver for the Euler equations
inline
void euler_rp(const real* q_left, const real* q_right,
              const real* aux_left, const real* aux_right,
              const void* aux_global_void,
              real* amdq, real* apdq, real* wave, real* s)
{
  const int num_states = euler_rp_grid_params.num_states;
  const euler_rp_aux_global_t* aux_global = (const euler_rp_aux_global_t *) aux_global_void;

  real a1,a2,a3,a4;
  real delta[4];
  real gamma1;
  real rhsqrtl,rhsqrtr,pl,pr,rhsq2,u,v,enth,u2v2,a,g1a2,euv;

  // Roe averages
  gamma1 = aux_global->gamma - 1.0;
  rhsqrtl = sqrt(q_right[0]);
  rhsqrtr = sqrt(q_left[0]);
  pl = gamma1*(q_right[3] - 0.5*(pow(q_right[1],2) + pow(q_right[2],2))/q_right[0]);
  pr = gamma1*(q_left[3] - 0.5*(pow(q_left[1],2) + pow(q_left[2],2))/q_left[0]);
  rhsq2 = rhsqrtl + rhsqrtr;
  u = (q_right[1]/rhsqrtl + q_left[1]/rhsqrtr) / rhsq2;
  v = (q_right[2]/rhsqrtl + q_left[2]/rhsqrtr) / rhsq2;
  enth = (((q_right[3]+pl)/rhsqrtl + (q_left[3]+pr)/rhsqrtr)) / rhsq2;
  u2v2 = pow(u,2) + pow(v,2);
  // This line has been changed to handle the random init conditions
  a2 = fabs(gamma1 * (enth - 0.5 * u2v2));
  a = sqrt(a2);
  g1a2 = gamma1 / a2;
  euv = enth - u2v2;

  // Wave strengths
  delta[0] = q_left[0] - q_right[0];
  delta[1] = q_left[1] - q_right[1];
  delta[2] = q_left[2] - q_right[2];
  delta[3] = q_left[3] - q_right[3];
  a3 = g1a2 * (euv*delta[0] + u*delta[1] + v*delta[2] - delta[3]);
  a2 = delta[2] - v*delta[0];
  a4 = (delta[1] + (a-u)*delta[0] - a*a3) / (2.0*a);
  a1 = delta[0] - a3 - a4;

  // Wave speeds
  s[0] = u-a;
  s[1] = u;
  s[2] = u;
  s[3] = u+a;

  // Waves
  wave[0*num_states + 0] = a1;
  wave[0*num_states + 1] = a1 * s[0];
  wave[0*num_states + 2] = a1 * v;
  wave[0*num_states + 3] = a1 * (enth - u * a);

  wave[1*num_states + 0] = 0.0;
  wave[1*num_states + 1] = 0.0;
  wave[1*num_states + 2] = a2;
  wave[1*num_states + 3] = a2 * v;

  wave[2*num_states + 0] = a3;
  wave[2*num_states + 1] = a3 * u;
  wave[2*num_states + 2] = a3 * v;
  wave[2*num_states + 3] = a3 * 0.5 * u2v2;

  wave[3*num_states + 0] = a4;
  wave[3*num_states + 1] = a4 * s[3];
  wave[3*num_states + 2] = a4 * v;
  wave[3*num_states + 3] = a4 * (enth + u * a);

  // Calculate fluctuations
  for (int eqn_index = 0; eqn_index < 4; eqn_index++)
  {
    amdq[eqn_index] = 0.0;
    apdq[eqn_index] = 0.0;
  }
  // Depending on the wave speed, amdq or apdq is modified
  for (int wave_index = 0; wave_index < 4; wave_index++)
  {
    if (s[wave_index] <= 0.0)
    {
      for (int eqn_index = 0; eqn_index < 4; eqn_index++)
      {
        amdq[eqn_index] += s[wave_index] * wave[wave_index + eqn_index*num_states];
      }
    }
    else
    {
      for (int eqn_index = 0; eqn_index < 4; eqn_index++)
      {
        apdq[eqn_index] += s[wave_index] * wave[wave_index + eqn_index*num_states];
      }
    }
  }
}

#endif
