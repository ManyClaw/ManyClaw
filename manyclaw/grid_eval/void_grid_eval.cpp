#include "void_grid_eval.h"

// Evaluates void_rp via serial execution
void void_rp_grid_eval_serial(rp_t rp,
                              rp_grid_params_t rp_grid_params,
                              const real* q,
                              const real* aux,
                              const void* aux_global,
                              const int nx,
                              const int ny,
                              real* amdq,
                              real* apdq,
                              real* wave,
                              real* wave_speed)
{
  unsigned col, row;
  const int num_ghost = rp_grid_params.num_ghost;
  const int num_eqn = rp_grid_params.num_eqn;
  const int num_aux = rp_grid_params.num_aux;
  const int num_wave = rp_grid_params.num_wave;

  FieldIndexer fi(nx, ny, num_ghost, num_eqn);
  FieldIndexer fi_aux(nx, ny, num_ghost, num_aux);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_wave);

  for(row = 1; row < efi.num_row_edge_transverse(); ++row){
    for(col = 1; col < efi.num_col_edge_normal() + 1; ++col) {
      rp(q + fi.idx(row, col-1), q + fi.idx(row, col),
         aux + fi_aux.idx(row, col-1), aux + fi_aux.idx(row, col), aux_global, 0,
         amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
         wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));

      rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
         aux + fi_aux.idx(row-1, col), aux + fi_aux.idx(row, col), aux_global, 1,
         amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
         wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
    }
  }

  for(col = 1; col < efi.num_col_edge_transverse(); ++col) {
    row = efi.num_row_edge_transverse();
    rp(q + fi.idx(row - 1, col), q + fi.idx(row, col),
       aux + fi_aux.idx(row - 1, col), aux + fi_aux.idx(row, col), aux_global, 1,
       amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
       wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
  }

  for(row = 1; row < efi.num_row_edge_transverse(); ++row){
    col = efi.num_col_edge_transverse();
    rp(q + fi.idx(row, col - 1), q + fi.idx(row, col),
       aux + fi_aux.idx(row, col - 1), aux + fi_aux.idx(row, col), aux_global, 0,
       amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
       wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));
  }
}

void void_rp_grid_eval_omp(rp_t rp,
                           rp_grid_params_t rp_grid_params,
                           const real* q,
                           const real* aux,
                           const void* aux_global,
                           const int nx,
                           const int ny,
                           real* amdq,
                           real* apdq,
                           real* wave,
                           real* wave_speed)
{
  unsigned col, row;
  const int num_ghost = rp_grid_params.num_ghost;
  const int num_eqn = rp_grid_params.num_eqn;
  const int num_aux = rp_grid_params.num_aux;
  const int num_wave = rp_grid_params.num_wave;

  FieldIndexer fi(nx, ny, num_ghost, num_eqn);
  FieldIndexer fi_aux(nx, ny, num_ghost, num_aux);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_wave);

#pragma omp parallel shared(q, aux, amdq, apdq, wave, wave_speed) private(col, row)
  {
#pragma omp for schedule(runtime) nowait
    for(row = 1; row < efi.num_row_edge_transverse(); ++row){
      for(col = 1; col < efi.num_col_edge_normal() + 1; ++col) {
        rp(q + fi.idx(row, col-1), q + fi.idx(row, col),
           aux + fi_aux.idx(row, col-1), aux + fi_aux.idx(row, col), aux_global, 0,
           amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
           wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));

        rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
           aux + fi_aux.idx(row-1, col), aux + fi_aux.idx(row, col), aux_global, 1,
           amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
           wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
      }
    }

#pragma omp for schedule(runtime) nowait
    for(col = 1; col < efi.num_col_edge_transverse(); ++col) {
      row = efi.num_row_edge_transverse();
      rp(q + fi.idx(row - 1, col), q + fi.idx(row, col),
         aux + fi_aux.idx(row - 1, col), aux + fi_aux.idx(row, col), aux_global, 1,
         amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
         wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
    }

#pragma omp for schedule(runtime) nowait
    for(row = 1; row < efi.num_row_edge_transverse(); ++row){
      col = efi.num_col_edge_transverse();
      rp(q + fi.idx(row, col - 1), q + fi.idx(row, col),
         aux + fi_aux.idx(row, col - 1), aux + fi_aux.idx(row, col), aux_global, 0,
         amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
         wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));
    }
  }
}

struct void_grid_eval_tbb_body
{
  rp_t rp;
  rp_grid_params_t rp_grid_params;    
  const real* q;
  const real* aux;
  const void* aux_global;
  const int nx;
  const int ny;
  real* amdq; real* apdq;
  real* wave; real* wave_speed;

  void_grid_eval_tbb_body(rp_t rp,
                          rp_grid_params_t rp_grid_params,
                          const real* q,  const real* aux, const void* aux_global,
                          const int nx, const int ny,
                          real* amdq, real* apdq, real* wave,
                          real* wave_speed)
    : rp(rp), rp_grid_params(rp_grid_params), q(q), aux(aux), aux_global(aux_global), nx(nx), ny(ny), amdq(amdq), apdq(apdq), wave(wave), wave_speed(wave_speed)
  {}

  void operator()(const tbb::blocked_range2d<int>& r) const
  {
    unsigned col, row;
    const int num_ghost = rp_grid_params.num_ghost;
    const int num_eqn = rp_grid_params.num_eqn;
    const int num_aux = rp_grid_params.num_aux;
    const int num_wave = rp_grid_params.num_wave;

    FieldIndexer fi(nx, ny, num_ghost, num_eqn);
    FieldIndexer fi_aux(nx, ny, num_ghost, num_aux);
    EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_wave);

    for(row = 1; row < efi.num_row_edge_transverse(); ++row){
      for(col = 1; col < efi.num_col_edge_normal() + 1; ++col) {
        rp(q + fi.idx(row, col-1), q + fi.idx(row, col),
           aux + fi_aux.idx(row, col-1), aux + fi_aux.idx(row, col), aux_global, 0,
           amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
           wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));

        rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
           aux + fi_aux.idx(row-1, col), aux + fi_aux.idx(row, col), aux_global, 1,
           amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
           wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
      }
    }

    for(col = 1; col < efi.num_col_edge_transverse(); ++col) {
      row = efi.num_row_edge_transverse();
      rp(q + fi.idx(row - 1, col), q + fi.idx(row, col),
         aux + fi_aux.idx(row - 1, col), aux + fi_aux.idx(row, col), aux_global, 1,
         amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
         wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
    }

    for(row = 1; row < efi.num_row_edge_transverse(); ++row){
      col = efi.num_col_edge_transverse();
      rp(q + fi.idx(row, col - 1), q + fi.idx(row, col),
         aux + fi_aux.idx(row, col - 1), aux + fi_aux.idx(row, col), aux_global, 0,
         amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
         wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));
    }
  }
};



void void_rp_grid_eval_tbb(rp_t rp,
                           rp_grid_params_t rp_grid_params,
                           const real* q,  const real* aux, const void* aux_global,
                           const int nx, const  int ny,
                           real* amdq, real* apdq, real* wave,
                           real* wave_speed)
{
  const int num_ghost = rp_grid_params.num_ghost;
  const int num_eqn = rp_grid_params.num_eqn;
  const int num_wave = rp_grid_params.num_wave;

  FieldIndexer fi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_wave);
  void_grid_eval_tbb_body body(rp, rp_grid_params, 
                               q, aux, aux_global, nx, ny, amdq, apdq, wave, wave_speed);

  // note: we use nx+1 and ny+1 here and < in the body (instead of <= in the serial reference)
  tbb::parallel_for(::tbb::blocked_range2d<int,int>(1, efi.num_row_edge_transverse(),
                                                    1, efi.num_col_edge_transverse()), body);
}
