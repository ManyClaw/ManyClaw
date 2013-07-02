#include "acoustics_const_rp_step_serial.h"
/**
Implements the serial stepping of the acoustics_const rieman problem through a structured
2D grid.

inputs: q, aux, numGhost, numStates, numWaves, nx, ny
outputs: amdq, apdq, wave, wave_speeds

 */
void acoustics_const_rp_step_serial(const real* q, const real* aux,
                              const int nx, const int ny, real* amdq, real* apdq,
                              real* wave, real* wave_speeds)
{
  int col, row, idx_left, idx_right, idx_up, idx_down, idx_out;
  const int num_ghost = acoustics_const_rp_grid_params.num_ghost;
  const int num_eqn = acoustics_const_rp_grid_params.num_eqn;
  const int num_waves = acoustics_const_rp_grid_params.num_waves;

  for(row = num_ghost; row <= ny + num_ghost; ++row)
  {
    for(col = num_ghost; col <= nx + num_ghost; ++col)
    {
      idx_left = col + row*(nx + 2*num_ghost) - 1;
      idx_right = idx_left + 1;
      idx_out = (col - num_ghost) + (row - num_ghost) * (nx + 1);
      acoustics_const_rp(q + idx_left*num_eqn, q + idx_right*num_eqn,
                   aux, aux, &acoustics_const_rp_aux_global,
                   amdq + idx_out*num_eqn, apdq + idx_out*num_eqn,
                   wave + num_waves*num_eqn*idx_out, wave_speeds + num_waves*idx_out);
    }
  }

  for(col = num_ghost; col <= nx + num_ghost; ++col)
  {
    for(row = num_ghost; row <= ny + num_ghost; ++row)
    {
      idx_up = col + (row - 1)*(nx + 2*num_ghost);
      idx_down = idx_up + nx + 2*num_ghost;
      idx_out = (col - num_ghost) + (row - num_ghost) * (nx + 1) + ((nx + 1)*(ny + 1));
      acoustics_const_rp(q + idx_up*num_eqn, q + idx_down*num_eqn,
                   aux, aux, &acoustics_const_rp_aux_global,
                   amdq + idx_out*num_eqn, apdq + idx_out*num_eqn,
                   wave + num_waves*num_eqn*idx_out, wave_speeds + num_waves*idx_out);
    }
  }
}

