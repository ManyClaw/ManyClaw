#include "advection_rp_step_serial_tiled.h"
#include <math.h>
#define TILE_SIZE 512

void advection_rp_step_serial_tiled( real* q,  real* aux,  int num_ghost,
   int num_states,  int num_waves,  int nx,  int ny,
  real* amdq, real* apdq, real* wave, real* wave_speeds)
{
  int col, row, idx_left, idx_center, idx_up, idx_out_x, idx_out_y;
  const int num_aux = 2;
  int row_tile, col_tile;

  for(row_tile=0; row_tile < ny/TILE_SIZE; ++row_tile){
    for(col_tile=0; col_tile < nx/TILE_SIZE; ++col_tile){
      for(row = row_tile * TILE_SIZE + num_ghost; row <= (row_tile+1)*TILE_SIZE + num_ghost; ++row) {
        for(col = col_tile * TILE_SIZE + num_ghost; col <= (col_tile+1)*TILE_SIZE + num_ghost; ++col) {
          idx_left = col + row*(nx + 2*num_ghost) - 1;
          idx_center = idx_left + 1;
          idx_out_x = (col - num_ghost) + (row - num_ghost) * (nx + 1);
          advection_rp(q + idx_left*num_states, q + idx_center*num_states, num_states,
                       aux + idx_left*num_aux, aux + idx_center*num_aux,
                       amdq + idx_out_x*num_states, apdq + idx_out_x*num_states,
                       wave + num_waves*num_states*idx_out_x, wave_speeds + num_waves*idx_out_x);
          idx_up = col + (row - 1)*(nx + 2*num_ghost);
          idx_out_y = (col - num_ghost) + (row - num_ghost) * (nx + 1) + (nx + 1)*(ny + 1);
          advection_rp(q + idx_up*num_states, q + idx_center*num_states, num_states,
                       aux + idx_up*num_aux, aux + idx_center*num_aux,
                       amdq + idx_out_y*num_states, apdq + idx_out_y*num_states,
                       wave + num_waves*num_states*idx_out_y, wave_speeds + num_waves*idx_out_y);
        }
      }
    }
  }
}

