#include "advection_rp_step_serial.h"
/**
Implements the serial stepping of the advection rieman problem through a structured
2D grid.

inputs: q, aux, numGhost, numStates, numWaves, nx, ny
outputs: amdq, apdq, wave, wave_speeds

 */
void advection_rp_step_serial( real* q,  real* aux,  int num_ghost,
   int num_states,  int num_waves,  int nx,  int ny,
  real* amdq, real* apdq, real* wave, real* wave_speeds)
{
  int col, row, idx_left, idx_right, idx_up, idx_down, idx_out;
  const int num_aux = 2;

  for(row = num_ghost; row <= ny + num_ghost; ++row)
  {
    for(col = num_ghost; col <= nx + num_ghost; ++col)
    {
      idx_left = col + row*(nx + 2*num_ghost) - 1;
      idx_right = idx_left + 1;
      idx_out = (col - num_ghost) + (row - num_ghost) * (nx + 1);
      advection_rp(q + idx_left*num_states, q + idx_right*num_states, num_states,
                   aux + idx_left*num_aux, aux + idx_right*num_aux,
                   amdq + idx_out*num_states, apdq + idx_out*num_states,
                   wave + num_waves*num_states*idx_out, wave_speeds + num_waves*idx_out);
    }
  }

  for(col = num_ghost; col <= nx + num_ghost; ++col)
  {
    for(row = num_ghost; row <= ny + num_ghost; ++row)
    {
      idx_up = col + (row - 1)*(nx + 2*num_ghost);
      idx_down = idx_up + nx + 2*num_ghost;
      idx_out = (col - num_ghost) + (row - num_ghost) * (nx + 1) + ((nx + 1)*(ny + 1));
      advection_rp(q + idx_up*num_states, q + idx_down*num_states, num_states,
                   aux + idx_up*num_aux, aux + idx_down*num_aux,
                   amdq + idx_out*num_states, apdq + idx_out*num_states,
                   wave + num_waves*num_states*idx_out, wave_speeds + num_waves*idx_out);
    }
  }
}

