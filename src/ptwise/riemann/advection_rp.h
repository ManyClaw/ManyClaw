#ifndef ADVECTION_RP_H
#define ADVECTION_RP_H

#include <math.h>
#include <algorithm>

typedef double real;

struct rp_params
{
  real u;
  int num_ghost;
  int num_states;
  int num_waves;
};

static const rp_params advection_rp_params = {1.0, 2, 1, 1};

// Constanst coefficient Riemann problem
inline
void advection_rp(const real* q_left, const real* q_right,
                  const real* aux_left, const real* aux_right,
                  const rp_params* params,
                  real* amdq, real* apdq, real* wave, real* s)
{
    // Wave and speed
    wave[0] = q_right[0] - q_left[0];
    s[0] = params->u;

    // This could also be implemented as an if statement or statically
    // based on the sign of s[0] (a.k.a. u)
    amdq[0] = std::min(0.0,s[0] * wave[0]);
    apdq[0] = std::max(0.0,s[0] * wave[0]);

}

typedef void (*advection_rp_step_t)(const real* q, const real* aux,
                                    const int nx, const int ny, real* amdq,  real* apdq,
                                    real* wave, real* wave_speeds);
#endif // ADVECTION_RP_H
