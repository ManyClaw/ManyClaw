#ifndef ADVECTION_RP_H
#define ADVECTION_RP_H

#include "../common/manyclaw_common.h"

struct advection_rp_aux_global_t
{
  real u[2];
};

static const advection_rp_aux_global_t advection_rp_aux_global_default = { {1.0, 1.0} };
static const rp_grid_params_t advection_rp_grid_params = {2, 1, 0, 1};

// Constanst coefficient Riemann problem
// TODO:  Deprecate advection_rp and change name to advection_const_rp
inline
void advection_rp(const real* q_left, const real* q_right,
                  const real* aux_left, const real* aux_right,
                  const void* aux_global_void,
                  const int direction,
                  real* amdq, real* apdq, real* wave, real* s)
{
  const advection_rp_aux_global_t* aux_global = (const advection_rp_aux_global_t *) aux_global_void;

  // Wave and speed
  wave[0] = q_right[0] - q_left[0];
  s[0] = aux_global->u[direction];

  // This could also be implemented as an if statement or statically
  // based on the sign of s[0] (a.k.a. u)
  amdq[0] = std::min(0.0, s[0]) * wave[0];
  apdq[0] = std::max(0.0, s[0]) * wave[0];
}

// Variable coefficient Riemann problem
static const rp_grid_params_t advection_var_rp_grid_params = {2, 1, 2, 1};

inline
void advection_var_rp(const real* q_left, const real* q_right,
                  const real* aux_left, const real* aux_right,
                  const void* aux_global_void,
                  const int direction,
                  real* amdq, real* apdq, real* wave, real* s)
{
  // Wave and speed
  wave[0] = q_right[0] - q_left[0];
  s[0] = 0.5 * (aux_left[direction] + aux_right[direction]);

  // This could also be implemented as an if statement or statically
  // based on the sign of s[0] (a.k.a. u)
  amdq[0] = std::min(0.0, s[0]) * wave[0];
  apdq[0] = std::max(0.0, s[0]) * wave[0];
}
#endif // ADVECTION_RP_H
