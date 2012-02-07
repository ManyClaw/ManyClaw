// Authors: Gorune Ohannessian
//          Andy R. Terrel

#include <math.h>
#include <stdlib.h>

typedef double real;

// Advection Riemann problem.
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

int advection_single_step(){
  // input states
  real q_left[3] = {0.0, 1.0, 1.0};
  real q_right[3] = {0.0, 1.25, 0.75};
  int numStates = 3;
  real u_left[2] = {1.0, 1.0};
  real u_right[2] = {1.2, 1.0};

  //output states
  real amdq[3];
  real apdq[3];
  real wave[3*3];
  real wave_speeds[3];
}

int advection_step_serial(real* q, real* aux, int nx, int ny,
                          int num_states) {
  real* waves = (real*) malloc(sizeof(real)*nx*ny*2*num_states);
  real* amdq = (real*) malloc(sizeof(real)*nx*ny*num_states);
  real* apdq = (real*) malloc(sizeof(real)*nx*ny*num_states);
  real* wave_speeds = (real*) malloc(sizeof(real)*nx*ny*2*num_states);
  int x_i, row, left, right;
  for(row = 0; row <= ny; ++row){
    for(x_i = 0; x_i <= nx; ++x_i){
      left = x_i + row * nx;
      right = left + 1;
      advection_rp(q + left, q + right, 3,
                   aux + left, aux + right,
                   amdq + left, apdq + left,
                   waves + 2*left, wave_speeds + 2*left);
    }
  }
  int y_i, col;
  for(col = 0; col <= nx; ++col){
    for(y_i = 0; y_i <= ny; ++y_i){
      left = y_i * nx + col;
      right = left + nx;
      advection_rp(q + left, q + right, 3,
                   aux + left, aux + right,
                   amdq + left, apdq + left,
                   waves + 2*left, wave_speeds + 2*left);
    }
  }
  free(waves);
}
