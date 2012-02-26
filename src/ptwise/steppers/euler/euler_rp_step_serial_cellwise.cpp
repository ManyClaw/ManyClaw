#include "euler_rp_step_serial_cellwise.h"

void euler_rp_step_serial_cellwise( const real* q,  const real* aux,
                                        const int nx, const  int ny,
                                        real* amdq, real* apdq, real* wave,
                                        real* wave_speeds)
{
  int col, row, idx_left, idx_center, idx_up, idx_out_x, idx_out_y;
  const int num_ghost = euler_rp_grid_params.num_ghost;
  const int num_states = euler_rp_grid_params.num_states;
  const int num_waves = euler_rp_grid_params.num_waves;

  for(row = num_ghost; row <= ny + num_ghost; ++row) {
    for(col = num_ghost; col <= nx + num_ghost; ++col) {
      idx_left = col + row*(nx + 2*num_ghost) - 1;
      idx_up = col + (row - 1)*(nx + 2*num_ghost);
      idx_center = idx_left + 1;
      idx_out_x = (col - num_ghost) + (row - num_ghost) * (nx + 1);
      idx_out_y = idx_out_x + ((nx + 1)*(ny + 1));
      euler_rp(q + idx_left*num_states, q + idx_center*num_states,
                   aux, aux,  &euler_rp_aux_global,
                   amdq + idx_out_x*num_states, apdq + idx_out_x*num_states,
                   wave + num_waves*num_states*idx_out_x, wave_speeds + num_waves*idx_out_x);
      euler_rp(q + idx_up*num_states, q + idx_center*num_states,
                   aux, aux,  &euler_rp_aux_global,
                   amdq + idx_out_y*num_states, apdq + idx_out_y*num_states,
                   wave + num_waves*num_states*idx_out_y, wave_speeds + num_waves*idx_out_y);
    }
  }
}
