#ifndef ADVECTION_RP_H
#define ADVECTION_RP_H

#include <math.h>

typedef double real;

// Advection Riemann problem.
inline
void advection_rp(real* q_left, real* q_right, int numStates,
                  real* u_left, real* u_right,	// input
                  real* amdq, real* apdq, real* wave,
                  real* wave_speeds) // output
{
  real rho_l = u_left[0];
  real bulk_l = u_left[1];

  real rho_r = u_right[0];
  real bulk_r = u_right[1];

  real c_l = sqrt(bulk_l/rho_l); // sound speed
  real z_l = c_l*rho_l; // impedance

  real c_r = sqrt(bulk_r/rho_r);
  real z_r = c_r*rho_r;

  // ( -(pr-pl) + zr(vr-vl) )/ (zl+zr)
  real alpha1 = ( q_left[0] - q_right[0] +
                   z_r*(q_right[1] - q_left[1])) / (z_l+z_r);
  // (  (pr-pl) + zl(vr-vl) )/ (zl+zr)
  real alpha2 = ( q_right[0] - q_left[0] +
                   z_l*(q_right[1] - q_left[1])) / (z_l+z_r);

  wave[0 + 0*numStates] = -alpha1*z_l;
  wave[1 + 0*numStates] = alpha1;
  wave[2 + 0*numStates] = 0;

  wave[0 + 1*numStates] = alpha2*z_r;
  wave[1 + 1*numStates] = alpha2;
  wave[2 + 1*numStates] = 0;

  amdq[0] = -c_l * wave[0 + 0*numStates];
  amdq[1] = -c_l * wave[1 + 0*numStates];
  amdq[2] =  0;	// 0   * wave[2 + 0*numStates];

  apdq[0] = c_r * wave[0 + 1*numStates];
  apdq[1] = c_r * wave[1 + 1*numStates];
  apdq[2] = 0;	// 0   * wave[2 + 1*numStates];
}

typedef void (*advection_rp_step_t)
  (real* q,     real* aux,
   int numGhost, int numStates, int numWaves, int nx, int ny, // inputs
   real* amdq,  real* apdq,
   real* wave,  real* wave_speeds);

#endif // ADVECTION_RP_H

