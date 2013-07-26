#ifdef TEMPLATE_GRID_EVAL_IPP
#define TEMPLATE_GRID_EVAL_IPP

// Macro this off since ForestClaw has trouble with templates.
#ifdef USE_TEMPLATE_GRID_EVAL
#include "template_grid_eval.h"

// Evaluates template_rp via serial execution
template<class rp_aux_global_t>
void template_rp_grid_eval_serial(rp_t rp,
                                  rp_grid_params_t rp_grid_params,
                                  const real* q,
                                  const real* aux,
                                  const rp_aux_global_t aux_global,
                                  const int nx,
                                  const int ny,
                                  real* amdq,
                                  real* apdq,
                                  real* wave,
                                  real* wave_speed)
{
  int col, row;
  const int num_ghost = rp_grid_params.num_ghost;
  const int num_eqn = rp_grid_params.num_eqn;
  const int num_wave = rp_grid_params.num_wave;

  FieldIndexer fi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_wave);

  for(row = 1; row < efi.num_row_edge_transverse(); ++row){
    for(col = 1; col < efi.num_col_edge_normal() + 1; ++col) {
      rp(q + fi.idx(row, col-1), q + fi.idx(row, col),
         aux, aux, aux_global, 0,
         amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
         wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));

      rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
         aux, aux, aux_global, 1,
         amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
         wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
    }
  }

  for(col = 1; col < efi.num_col_edge_transverse(); ++col) {
    row = efi.num_row_edge_transverse();
    rp(q + fi.idx(row - 1, col), q + fi.idx(row, col),
       aux, aux, aux_global, 1,
       amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
       wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
  }

  for(row = 1; row < efi.num_row_edge_transverse(); ++row){
    col = efi.num_col_edge_transverse();
    rp(q + fi.idx(row, col - 1), q + fi.idx(row, col),
       aux, aux, aux_global, 0,
       amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
       wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));
  }
}

#endif // USE_TEMPLATE_GRID_EVAL
#endif // TEMPLATE_GRID_EVAL_IPP

