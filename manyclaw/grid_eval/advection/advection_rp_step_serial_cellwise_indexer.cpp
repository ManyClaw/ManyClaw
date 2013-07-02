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

  for(row = 1; row < efi.num_row_edge_transverse(); ++row){
    for(col = 1; col < efi.num_col_edge_normal() + 1; ++col) {
      // printf("-> Calling %d, %d ", row, col);
      // printf("-> fi.idx %d, efi left %d down %d\n", fi.idx(row, col), efi.left_edge(row, col), efi.down_edge(row, col));
      advection_rp(q + fi.idx(row, col-1), q + fi.idx(row, col),
                   aux, aux,  &advection_rp_aux_global,
                   amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
                   wave + efi.left_edge(row, col), wave_speeds + efi.left_edge(row, col));

        advection_rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
                   aux, aux,  &advection_rp_aux_global,
                   amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
                     wave + efi.down_edge(row, col), wave_speeds + efi.down_edge(row, col));
    }
  }

  for(col = 1; col < efi.num_col_edge_transverse(); ++col) {
    row = efi.num_row_edge_transverse();
    // printf("-> Calling %d, %d ", row, col);
    // printf("-> fi.idx %d, efi left %d down %d\n", fi.idx(row, col), efi.left_edge(row, col), efi.down_edge(row, col));
    advection_rp(q + fi.idx(row - 1, col), q + fi.idx(row, col),
                 aux, aux,  &advection_rp_aux_global,
                 amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
                 wave + efi.down_edge(row, col), wave_speeds + efi.down_edge(row, col));
  }

  for(row = 1; row < efi.num_row_edge_transverse(); ++row){
    col = efi.num_col_edge_transverse();
    // printf("-> Calling %d, %d ", row, col);
    // printf("-> fi.idx %d, efi left %d down %d\n", fi.idx(row, col), efi.left_edge(row, col), efi.down_edge(row, col));
    advection_rp(q + fi.idx(row, col - 1), q + fi.idx(row, col),
                 aux, aux,  &advection_rp_aux_global,
                 amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
                 wave + efi.left_edge(row, col), wave_speeds + efi.left_edge(row, col));
  }
}

