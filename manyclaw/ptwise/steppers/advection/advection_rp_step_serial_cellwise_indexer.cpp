#include "advection_rp_step_serial_cellwise_indexer.h"

void advection_rp_step_serial_cellwise_indexer( const real* q,  const real* aux,
                                        const int nx, const  int ny,
                                        real* amdq, real* apdq, real* wave,
                                        real* wave_speeds)
{
  int col, row;
  const int num_ghost = advection_rp_grid_params.num_ghost;
  const int num_eqn = advection_rp_grid_params.num_eqn;
  const int num_waves = advection_rp_grid_params.num_waves;

  FieldIndexer fi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_waves);

  for(row = 0; row < nx + 2*num_ghost - 1; ++row){
    for(col = 0; col <= nx + 2*num_ghost - 1; ++col) {
      advection_rp(q + fi.idx(row, col), q + fi.idx(row, col + 1),
                   aux, aux,  &advection_rp_aux_global,
                   amdq + efi.right_edge(row, col), apdq + efi.right_edge(row, col),
                   wave + efi.right_edge(row, col), wave_speeds + efi.right_edge(row, col));

        advection_rp(q + fi.idx(row, col), q + fi.idx(row + 1, col),
                   aux, aux,  &advection_rp_aux_global,
                   amdq + efi.up_edge(row, col), apdq + efi.up_edge(row, col),
                     wave + efi.up_edge(row, col), wave_speeds + efi.up_edge(row, col));
    }
  }
}

