#ifndef ACOUSTICS_VAR_RP_H
#define ACOUSTICS_VAR_RP_H

#include "../../common/manyclaw_common.h"

struct acoustics_var_rp_aux_global_t
{
};

static const acoustics_var_rp_aux_global_t acoustics_var_rp_aux_global = {};
static const rp_grid_params acoustics_var_rp_grid_params = {2, 3, 2, 2};

// Point-wise constant acoustics Riemann solver
// aux_* = [rho,bulk]
inline
void acoustics_var_rp(const real* q_left, const real* q_right,
                      const real* aux_left, const real* aux_right,
                      const acoustics_var_rp_aux_global_t* aux_global,
                      real* amdq, real* apdq, real *wave, real *s)
{
    const int num_eqn = acoustics_var_rp_grid_params.num_eqn;

    // Local physical constants
    real c_left = sqrt(aux_left[1] / aux_left[0]);
    real c_right = sqrt(aux_right[1] / aux_right[0]);
    real Z_left = c_left * aux_left[0];
    real Z_right = c_right * aux_right[0];

    // Wave strengths
    real delta[2],alpha[2];
    delta[0] = q_right[0] - q_left[0];
    delta[1] = q_left[1] - q_right[1];
    alpha[0] = (-delta[0] + Z_right * delta[1]) / (Z_right + Z_left);
    alpha[1] = ( delta[0] + Z_left * delta[1])  / (Z_right + Z_left);

    // Wave speeds
    s[0] = -c_left;
    s[1] =  c_right;

    // Set wave vectors
    wave[0 * num_eqn + 0] = -alpha[0] * Z_left;
    wave[0 * num_eqn + 1] = alpha[0];
    wave[0 * num_eqn + 2] = 0.0;

    wave[1 * num_eqn + 0] = alpha[1] * Z_right;
    wave[1 * num_eqn + 1] = alpha[1];
    wave[1 * num_eqn + 2] = 0.0;

    // Grid edge fluctuations
    amdq[0] = s[0] * wave[0 * num_eqn + 0];
    amdq[1] = s[0] * wave[0 * num_eqn + 1];
    amdq[2] = s[0] * wave[0 * num_eqn + 2];  // This could just be set to zero

    apdq[0] = s[1] * wave[1 * num_eqn + 0];
    apdq[1] = s[1] * wave[1 * num_eqn + 1];
    apdq[2] = s[1] * wave[1 * num_eqn + 2];  // This could just be set to zero
}

#endif // ACOUSTICS_VAR_RP_H
