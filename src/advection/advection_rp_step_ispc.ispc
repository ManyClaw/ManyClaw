typedef double real;

// Advection Riemann problem.
void advection_rp(real* q_left, real* q_right, uniform int numStates,
                  real* u_left, real* u_right,	// input
                  real* amdq, real* apdq,
                  real* wave, real* wave_speeds) // output
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

export void advection_rp_step_ispc
  (uniform real* q,     uniform real* aux,
   uniform int numGhost, uniform int numStates, uniform int numWaves, uniform int nx, uniform int ny, // inputs
   uniform real* amdq,  uniform real* apdq,
   uniform real* wave,  uniform real* wave_speeds)
{
  foreach(row = 0 ... ny, x_i = 0 ... nx)
  {
      int left = x_i + row * nx;
      int right = left + 1;
      advection_rp(q + left, q + right, numStates,
                   aux + left, aux + right,
                   amdq + left, apdq + left,
                   wave + 2*left, wave_speeds + 2*left);
  }

  foreach(col = 0 ... nx, y_i = 0 ... ny)
  {
      int left = y_i * nx + col;
      int right = left + nx;
      advection_rp(q + left, q + right, numStates,
                   aux + left, aux + right,
                   amdq + left, apdq + left,
                   wave + 2*left, wave_speeds + 2*left);
  }
}

//export void advection_step_ispc_task
//        (uniform int N,
//         uniform real* q_left, uniform real* q_right, uniform int numStates,
//         uniform real* u_left, uniform real* u_right,	// input
//         uniform real* amdq,   uniform real* apdq, 
//         uniform real* wave,
//         uniform real* wave_speeds) // output
//{
//    advection_step_ispc
//        (N,
//         q_left + taskIndex, q_right + taskIndex, numStates + taskIndex,
//         u_left + taskIndex, u_right + taskIndex,
//         amdq + taskIndex, apdq + taskIndex,
//         wave + taskIndex, wave_speeds + taskIndex);
//}
//
//export void advection_step_ispc_tasks
//        (uniform int N,
//         uniform real* q_left, uniform real* q_right, uniform int numStates,
//         uniform real* u_left, uniform real* u_right,	// input
//         uniform real* amdq,   uniform real* apdq, 
//         uniform real* wave,
//         uniform real* wave_speeds) // output
//{
//    launch[N] < advection_step_ispc
//        (N,
//         q_left, q_right, numStates,
//         u_left, u_right,
//         amdq, apdq,
//         wave, wave_speeds) >;
//
//     sync;
//}
//
