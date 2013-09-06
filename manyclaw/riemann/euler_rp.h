#ifndef EULER_RP_H
#define EULER_RP_H

#include "../common/manyclaw_common.h"

#include <iostream>

struct euler_rp_aux_global_t
{
    real gamma;
    bool entropy_fix;
};

static const euler_rp_aux_global_t euler_rp_aux_global_default = {1.4, false};
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

  // Apply entropy fix
  if (aux_global->entropy_fix)
  {
    // TODO: Remove statement, only here for debugging
    if (delta[0] != 0.0 && false)
    {
      std::cout << "Non trivial problem\n";
      std::cout << " q_left = [" << q_left[0] << ", " << q_left[1] << ", " << q_left[2] << ", " << q_left[3] << "]\n";
      std::cout << " q_right = [" << q_right[0] << ", " << q_right[1] << ", " << q_right[2] << ", " << q_right[3] << "]\n";
      std::cout << " wave[0] = [" << wave[0*num_eqn + 0] << ", " << wave[0*num_eqn + 1] << ", " << wave[0*num_eqn + 2] << ", " << wave[0*num_eqn + 3] << "]\n";
      std::cout << " wave[3] = [" << wave[3*num_eqn + 0] << ", " << wave[3*num_eqn + 1] << ", " << wave[3*num_eqn + 2] << ", " << wave[3*num_eqn + 3] << "]\n";
    }

    // Splitting storage
    real sfract, df;

    // Left state (original input to left of cell interface)
    real rhoim1, pim1, cim1, s0;

    // Right state (original input to right of cell interface)
    real rhoi, pi, ci, s3;

    // State to the right of 1-wave
    real rho1, rhou1, rhov1, E1, p1, c1, s1; 

    // State to the left of 4-wave
    real rho2, rhou2, rhov2, E2, p2, c2, s2;

    // Check 1-wave
    rhoim1 = q_left[0];
    pim1 = (gamma - 1.0) * (q_left[3] - 0.5 * (pow(q_left[normal_velocity], 2) + pow(q_left[trans_velocity], 2)) / rhoim1);
    cim1 = sqrt(gamma * pim1 / rhoim1);
    s0 = q_left[normal_velocity] / rhoim1 - cim1; // u-c in left state

    // Check for fully supersonic case
    if (s0 >= 0.0 && s[0] > 0.0)
    {
      // Everything is right-going
      for (int eqn_index = 0; eqn_index < num_eqn; ++eqn_index)
      {
        amdq[eqn_index] = 0.0;
      }
    }
    else
    {
      // Check what speed is on the left side of the 1-wave
      rho1 = q_left[0] + wave[0 * num_eqn + 0];
      rhou1 = q_left[normal_velocity] + wave[0 * num_eqn + normal_velocity];
      rhov1 = q_left[trans_velocity] + wave[0 * num_eqn + trans_velocity];
      E1 = q_left[3] + wave[0 * num_eqn + 3];
      p1 = (gamma - 1.0) * (E1 - 0.5 * (pow(rhou1, 2) + pow(rhov1, 2)) / rho1);
      c1 = sqrt(gamma * p1 / rho1);
      s1 = rhou1 / rho1 - c1; // u-c to right of 1-wave
      if (s0 < 0.0 && s1 > 0.0)
      {
        // Transonic rarefaction in the 1-wave
        sfract = s0 * (s1 - s[0]) / (s1 - s0);
      }
      else if (s[0] < 0.0)
      {
        // 1-wave is left-going
        sfract = s[0];
      }
      else
      {
        // 1-wave is right-going, this should not happen since s0 < 0
        sfract = 0.0;  
      }

      // Accumulate fractions of fluctuations that are left going
      for (int eqn_index = 0; eqn_index < num_eqn; ++eqn_index)
      {
        amdq[eqn_index] = sfract * wave[0 * num_eqn + eqn_index];
      }

      // Check contact discontinuity
      // If this is not true than 2- and 3-waves are rightgoing
      if (s[1] < 0.0) 
      {
        for (int eqn_index = 0; eqn_index < num_eqn; ++eqn_index)
        {
          amdq[eqn_index] += s[1] * wave[1 * num_eqn + eqn_index];
          apdq[eqn_index] += s[2] * wave[2 * num_eqn + eqn_index];
        }

        // Check 4-wave
        rhoi = q_right[0];
        pi = (gamma - 1.0) * (q_right[3] - 0.5 * (pow(q_right[normal_velocity], 2) + pow(q_right[trans_velocity], 2)) / rhoi);
        ci = sqrt(gamma * pi / rhoi);
        s3 = q_right[normal_velocity] / rhoi + ci; // u+c in right state

        rho2 = q_right[0] - wave[3 * num_eqn + 0];
        rhou2 = q_right[normal_velocity] - wave[3 * num_eqn + normal_velocity];
        rhov2 = q_right[trans_velocity] - wave[3 * num_eqn + trans_velocity];
        E2 = q_right[3] - wave[3 * num_eqn + 3];
        p2 = (gamma - 1.0) * (E2 - 0.5 * (pow(rhou2, 2) + pow(rhov2, 2)) / rho2);
        c2 = sqrt(gamma * p2 / rho2);
        s2 = rhou2 / rho2 + c2; // u+c to left of 4-wave
        if (s2 < 0.0 && s3 > 0.0)
        {
          // Transonic rarefaction in the 4-wave
          sfract = s2 * (s3 - s[3]) / (s3 - s2);
        }
        else if (s[3] < 0.0)
        {
          // 4-wave is left going
          sfract = s[3];
        }
        else
        {
          // 4-wave is right going
          // Do not really need to do this but this is easier to do
          sfract = 0.0;
        }

        for (int eqn_index = 0; eqn_index < num_eqn; ++eqn_index)
        {
          amdq[eqn_index] += sfract * wave[3 * num_eqn + eqn_index];
        }
      }
    }

    // Compute the right going flux differences:
    // df = SUM s * wave    is the total flux difference and apdq = df - amdq
    for (int eqn_index = 0; eqn_index < num_eqn; ++eqn_index)
    {
      df = 0.0;
      for (int wave_index = 0; wave_index < num_wave; ++wave_index)
      {
        df += s[wave_index] * wave[wave_index * num_eqn + eqn_index];
      }
      apdq[eqn_index] = df - amdq[eqn_index];
    }

  }
  // No entropy fix
  else
  {
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


}

#endif
