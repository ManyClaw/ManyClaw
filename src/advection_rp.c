// Authors: Gorune Ohannessian
//          Andy R. Terrel

#include <math.h>
typedef double real;

// Advection Riemann problem.
void rp_advecion(real* q_left, real* q_right, int numStates,
                 real* u_left, real* u_right,	// input
                 real* amdq, real* apdq, real* wave,
                 real* wave_speeds) // output
{
  float rho_l = u_left[0];
  float bulk_l = u_left[1];

  float rho_r = u_right[0];
  float bulk_r = u_right[1];

  float c_l = sqrt(bulk_l/rho_l); // sound speed
  float z_l = c_l*rho_l; // impedance

  float c_r = sqrt(bulk_r/rho_r);
  float z_r = c_r*rho_r;

  // ( -(pr-pl) + zr(vr-vl) )/ (zl+zr)
  float alpha1 = ( q_left[0] - q_right[0] +
                   z_r*(q_right[1] - q_left[1])) / (z_l+z_r);
  // (  (pr-pl) + zl(vr-vl) )/ (zl+zr)
  float alpha2 = ( q_right[0] - q_left[0] +
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


int main(int argc, char** argv) {
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
