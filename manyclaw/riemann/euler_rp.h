#ifndef EULER_RP_H
#define EULER_RP_H

#include "../common/manyclaw_common.h"

struct euler_rp_aux_global_t
{
    real gamma;
};

static const euler_rp_aux_global_t euler_rp_aux_global_default = {1.4};
static const rp_grid_params_t euler_rp_grid_params = {2, 4, 0, 4};

// Implementation of the Roe solver for the Euler equations
inline
void euler_rp(const real* q_left, const real* q_right,
              const real* aux_left, const real* aux_right,
              const void* aux_global_void,
              const int direction,
              real* amdq, real* apdq, real* wave, real* s)
{
  const int num_eqn = euler_rp_grid_params.num_eqn;
  const int num_wave = euler_rp_grid_params.num_wave;
  const euler_rp_aux_global_t* aux_global = (const euler_rp_aux_global_t *) aux_global_void;


  // Direction of sweep
  int normal_velocity, trans_velocity;
  if (direction == 0)
  {
    normal_velocity = 1;
    trans_velocity = 2;
  }
  else
  {
    normal_velocity = 2;
    trans_velocity = 1;
  }

  // Compute Roe averages
  real gamma, p_l, p_r, u_hat, v_hat, H;
  gamma = aux_global->gamma;

  p_l = (gamma - 1.0) * (q_left[3] 
                   - 0.5 * (pow(q_left[1], 2) + pow(q_left[2], 2)) / q_left[0]);
  p_r = (gamma - 1.0) * (q_right[3] 
                - 0.5 * (pow(q_right[1], 2) + pow(q_right[2], 2)) / q_right[0]);

  u_hat = (q_left[normal_velocity] / sqrt(q_left[0]) 
                        + q_right[normal_velocity] / sqrt(q_right[0])) 
                                         / (sqrt(q_left[0]) + sqrt(q_right[0]));
  v_hat = (q_left[trans_velocity] / sqrt(q_left[0]) 
                        + q_right[trans_velocity] / sqrt(q_right[0])) 
                                         / (sqrt(q_left[0]) + sqrt(q_right[0]));
  H = ((q_left[3] + p_l) / sqrt(q_left[0]) 
                        + (q_right[3] + p_r) / sqrt(q_right[0])) 
                                         / (sqrt(q_left[0]) + sqrt(q_right[0]));
  
  // Wave strengths
  real alpha[4];
  real delta[4];
  delta[0] = q_right[0] - q_left[0];
  delta[1] = q_right[normal_velocity] - q_left[normal_velocity];
  delta[2] = q_right[trans_velocity] - q_left[trans_velocity];
  delta[3] = q_right[3] - q_left[3];

  real u2v2, c2;
  u2v2 = pow(u_hat, 2) + pow(v_hat, 2);
  c2 = (gamma - 1.0) * (H - 0.5 * u2v2);

  alpha[2] = (gamma - 1.0) / c2 * ((H - u2v2) * delta[0] + u_hat * delta[1] + v_hat * delta[2] - delta[3]);
  alpha[1] = delta[2] - v_hat * delta[0];
  alpha[3] = (delta[1] + (sqrt(c2) - u_hat) * delta[0] - sqrt(c2) * alpha[2]) / (2.0 * sqrt(c2));
  alpha[0] = delta[0] - alpha[2] - alpha[3];

  // Wave speeds
  s[0] = u_hat - sqrt(c2);
  s[1] = u_hat;
  s[2] = u_hat;
  s[3] = u_hat + sqrt(c2);

  // Waves
  wave[0 * num_eqn + 0] = alpha[0];
  wave[0 * num_eqn + 1] = alpha[0] * s[0];
  wave[0 * num_eqn + 2] = alpha[0] * v_hat;
  wave[0 * num_eqn + 3] = alpha[0] * (H - u_hat * sqrt(c2));

  wave[1 * num_eqn + 0] = 0.0;
  wave[1 * num_eqn + 1] = 0.0;
  wave[1 * num_eqn + 2] = alpha[1];
  wave[1 * num_eqn + 3] = alpha[1] * v_hat;

  wave[2 * num_eqn + 0] = alpha[2];
  wave[2 * num_eqn + 1] = alpha[2] * u_hat;
  wave[2 * num_eqn + 2] = alpha[2] * v_hat;
  wave[2 * num_eqn + 3] = alpha[2] * 0.5 * u2v2;

  wave[3 * num_eqn + 0] = alpha[3];
  wave[3 * num_eqn + 1] = alpha[3] * s[3];
  wave[3 * num_eqn + 2] = alpha[3] * v_hat;
  wave[3 * num_eqn + 3] = alpha[3] * (H + u_hat * sqrt(c2));

  // Calculate fluctuations
  for (int eqn_index = 0; eqn_index < num_eqn; ++eqn_index)
  {
    amdq[eqn_index] = 0.0;
    apdq[eqn_index] = 0.0;
  }

  for (int wave_index = 0; wave_index < num_wave; ++wave_index)
  {
    if (s[wave_index] <= 0.0)
    {
      for (int eqn_index = 0; eqn_index < num_eqn; ++eqn_index)
        amdq[eqn_index] += s[wave_index] * wave[eqn_index + wave_index * num_eqn];
    }
    else
    {
      for (int eqn_index = 0; eqn_index < num_eqn; ++eqn_index)
        apdq[eqn_index] += s[wave_index] * wave[eqn_index + wave_index * num_eqn];
    }
  }
}

#endif
