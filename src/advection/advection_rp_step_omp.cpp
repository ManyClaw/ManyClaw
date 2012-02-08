#include "advection_rp_step_omp.h"


void advection_rp_step_omp(
real* q,     real* aux,
int numGhost, int numStates, int numWaves, int nx, int ny, // inputs
real* amdq,  real* apdq,
real* wave,  real* wave_speeds
)
{
  int x_i, row, left, right, y_i, col;

#pragma omp parallel shared(q, aux, amdq, apdq, wave, wave_speeds) private(x_i, row, left, right, y_i, col)
  {
#pragma omp for schedule(runtime) nowait
    for(row = 0; row <= ny; ++row)
    {
      for(x_i = 0; x_i <= nx; ++x_i)
      {
        left = x_i + row * nx;
        right = left + 1;
        advection_rp(q + left, q + right, 3,
                     aux + left, aux + right,
                     amdq + left, apdq + left,
                     wave + 2*left, wave_speeds + 2*left);
      }
    }

#pragma omp for schedule(runtime) nowait
    for(col = 0; col <= nx; ++col)
    {
      for(y_i = 0; y_i <= ny; ++y_i)
      {
        left = y_i * nx + col;
        right = left + nx;
        advection_rp(q + left, q + right, 3,
                     aux + left, aux + right,
                     amdq + left, apdq + left,
                     wave + 2*left, wave_speeds + 2*left);
      }
    }
  }
}
