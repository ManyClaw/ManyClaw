#include "void_grid_eval.h"

// TODO: Remove cout
#include <iostream>

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
  const unsigned num_ghost = rp_grid_params.num_ghost;
  const unsigned num_eqn = rp_grid_params.num_eqn;
  const unsigned num_aux = rp_grid_params.num_aux;
  const unsigned num_wave = rp_grid_params.num_wave;

  FieldIndexer fi(nx, ny, num_ghost, num_eqn);
  FieldIndexer fi_aux(nx, ny, num_ghost, num_aux);

  EdgeFieldIndexer efi_wave(nx, ny, num_ghost, num_eqn, num_wave);
  EdgeFieldIndexer efi_asdq(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer efi_s(nx, ny, num_ghost, num_wave);

  // Interior of domain minus the top and right boundaries
  for(row = num_ghost; row < ny + num_ghost; ++row) {
    for(col = num_ghost; col < nx + num_ghost; ++col) {

      // Left edge of cell
      rp(q + fi.idx(row, col-1), q + fi.idx(row, col),
         aux + fi_aux.idx(row, col-1), aux + fi_aux.idx(row, col), aux_global, 0,
         amdq + efi_asdq.left_edge(row, col), apdq + efi_asdq.left_edge(row, col),
         wave + efi_wave.left_edge(row, col), wave_speed + efi_s.left_edge(row, col));

      // Bottom edge of cell
      rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
         aux + fi_aux.idx(row-1, col), aux + fi_aux.idx(row, col), aux_global, 1,
         amdq + efi_asdq.down_edge(row, col), apdq + efi_asdq.down_edge(row, col),
         wave + efi_wave.down_edge(row, col), wave_speed + efi_s.down_edge(row, col));
    }
  }

  // Bottom boundary ghost cells
  for(row = 0; row < num_ghost - 1; ++row) {
    for(col = num_ghost; col < nx + num_ghost; ++col) {

      rp(q + fi.idx(row, col), q + fi.idx(row+1, col),
         aux + fi_aux.idx(row-1, col), aux + fi_aux.idx(row, col), aux_global, 1,
         amdq + efi_asdq.up_edge(row, col), apdq + efi_asdq.up_edge(row, col),
         wave + efi_wave.up_edge(row, col), wave_speed + efi_s.up_edge(row, col));
    }  
  }

  // Left boundary ghost cells
  for(row = num_ghost; row < ny + num_ghost; ++row) {
    for (col = 0; col < num_ghost - 1; ++col) {

      rp(q + fi.idx(row, col), q + fi.idx(row+1, col),
         aux + fi_aux.idx(row, col), aux + fi_aux.idx(row, col + 1), aux_global, 1,
         amdq + efi_asdq.right_edge(row, col), apdq + efi_asdq.right_edge(row, col),
         wave + efi_wave.right_edge(row, col), wave_speed + efi_s.right_edge(row, col));
    }
  }

  // Right boundary ghost cells + right boundary edge
  for(row = num_ghost; row < ny + num_ghost; ++row) {
    for (col = nx + num_ghost; col < nx + 2 * num_ghost; ++col) {

      rp(q + fi.left(row, col), q + fi.idx(row, col),
         aux + fi_aux.idx(row, col), aux + fi_aux.idx(row, col + 1), aux_global, 1,
         amdq + efi_asdq.left_edge(row, col), apdq + efi_asdq.left_edge(row, col),
         wave + efi_wave.left_edge(row, col), wave_speed + efi_s.left_edge(row, col));
    }
  }

  // Top boundary ghost cells + top boundary edge
  for(row = ny + num_ghost; row < ny + 2 * num_ghost; ++row) {
    for (col = num_ghost; col < nx + num_ghost; ++col) {

      rp(q + fi.idx(row - 1, col), q + fi.idx(row - 1, col),
         aux + fi_aux.idx(row, col), aux + fi_aux.idx(row, col + 1), aux_global, 1,
         amdq + efi_asdq.down_edge(row, col), apdq + efi_asdq.down_edge(row, col),
         wave + efi_wave.down_edge(row, col), wave_speed + efi_s.down_edge(row, col));
    }
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
  
  EdgeFieldIndexer efi_wave(nx, ny, num_ghost, num_eqn, num_wave);
  EdgeFieldIndexer efi_asdq(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer efi_s(nx, ny, num_ghost, num_wave);

#pragma omp parallel shared(q, aux, amdq, apdq, wave, wave_speed) private(col, row)
  {
#pragma omp for schedule(runtime) nowait
    for(row = 1; row < efi_wave.num_row_edge_transverse(); ++row){
      for(col = 1; col < efi_wave.num_col_edge_normal() + 1; ++col) {
        rp(q + fi.idx(row, col-1), q + fi.idx(row, col),
           aux + fi_aux.idx(row, col-1), aux + fi_aux.idx(row, col), aux_global, 0,
           amdq + efi_asdq.left_edge(row, col), apdq + efi_asdq.left_edge(row, col),
           wave + efi_wave.left_edge(row, col), wave_speed + efi_s.left_edge(row, col));

        rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
           aux + fi_aux.idx(row-1, col), aux + fi_aux.idx(row, col), aux_global, 1,
           amdq + efi_asdq.down_edge(row, col), apdq + efi_asdq.down_edge(row, col),
           wave + efi_wave.down_edge(row, col), wave_speed + efi_s.down_edge(row, col));
      }
    }

#pragma omp for schedule(runtime) nowait
    for(col = 1; col < efi_wave.num_col_edge_transverse(); ++col) {
      row = efi_wave.num_row_edge_transverse();
      rp(q + fi.idx(row - 1, col), q + fi.idx(row, col),
         aux + fi_aux.idx(row - 1, col), aux + fi_aux.idx(row, col), aux_global, 1,
         amdq + efi_asdq.down_edge(row, col), apdq + efi_asdq.down_edge(row, col),
         wave + efi_wave.down_edge(row, col), wave_speed + efi_s.down_edge(row, col));
    }

#pragma omp for schedule(runtime) nowait
    for(row = 1; row < efi_wave.num_row_edge_transverse(); ++row){
      col = efi_wave.num_col_edge_transverse();
      rp(q + fi.idx(row, col - 1), q + fi.idx(row, col),
         aux + fi_aux.idx(row, col - 1), aux + fi_aux.idx(row, col), aux_global, 0,
         amdq + efi_asdq.left_edge(row, col), apdq + efi_asdq.left_edge(row, col),
         wave + efi_wave.left_edge(row, col), wave_speed + efi_s.left_edge(row, col));
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

  void operator()(const tbb::blocked_range2d<unsigned>& r) const
  {
    unsigned col, row;
    const int num_ghost = rp_grid_params.num_ghost;
    const int num_eqn = rp_grid_params.num_eqn;
    const int num_aux = rp_grid_params.num_aux;
    const int num_wave = rp_grid_params.num_wave;

    FieldIndexer fi(nx, ny, num_ghost, num_eqn);
    FieldIndexer fi_aux(nx, ny, num_ghost, num_aux);
  
    EdgeFieldIndexer efi_wave(nx, ny, num_ghost, num_eqn, num_wave);
    EdgeFieldIndexer efi_asdq(nx, ny, num_ghost, num_eqn);
    EdgeFieldIndexer efi_s(nx, ny, num_ghost, num_wave);

    // Should there be an efi_wave.num_row_edge_normal() here?
    for(row = r.rows().begin(); row != r.rows().end(); ++row){
      for(col = r.cols().begin(); col != std::min(r.cols().end(), efi_wave.num_col_edge_normal() + 1); ++col) {
        rp(q + fi.idx(row, col-1), q + fi.idx(row, col),
           aux + fi_aux.idx(row, col-1), aux + fi_aux.idx(row, col), aux_global, 0,
           amdq + efi_asdq.left_edge(row, col), apdq + efi_asdq.left_edge(row, col),
           wave + efi_wave.left_edge(row, col), wave_speed + efi_s.left_edge(row, col));

        rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
           aux + fi_aux.idx(row-1, col), aux + fi_aux.idx(row, col), aux_global, 1,
           amdq + efi_asdq.down_edge(row, col), apdq + efi_asdq.down_edge(row, col),
           wave + efi_wave.down_edge(row, col), wave_speed + efi_s.down_edge(row, col));
      }
    }
  }
};


struct void_grid_eval_tbb_boundary
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

  void_grid_eval_tbb_boundary(rp_t rp,
                          rp_grid_params_t rp_grid_params,
                          const real* q,  const real* aux, const void* aux_global,
                          const int nx, const int ny,
                          real* amdq, real* apdq, real* wave,
                          real* wave_speed)
    : rp(rp), rp_grid_params(rp_grid_params), q(q), aux(aux), aux_global(aux_global), nx(nx), ny(ny), amdq(amdq), apdq(apdq), wave(wave), wave_speed(wave_speed)
  {}

  void operator()(const tbb::blocked_range<unsigned>& r) const
  {
    unsigned col, row;
    const int num_ghost = rp_grid_params.num_ghost;
    const int num_eqn = rp_grid_params.num_eqn;
    const int num_aux = rp_grid_params.num_aux;
    const int num_wave = rp_grid_params.num_wave;

    FieldIndexer fi(nx, ny, num_ghost, num_eqn);
    FieldIndexer fi_aux(nx, ny, num_ghost, num_aux);
  
    EdgeFieldIndexer efi_wave(nx, ny, num_ghost, num_eqn, num_wave);
    EdgeFieldIndexer efi_asdq(nx, ny, num_ghost, num_eqn);
    EdgeFieldIndexer efi_s(nx, ny, num_ghost, num_wave);

    for(col = r.begin(); col != std::min(r.end(), efi_wave.num_col_edge_transverse()); ++col) {
      row = efi_wave.num_row_edge_transverse();
      rp(q + fi.idx(row - 1, col), q + fi.idx(row, col),
         aux + fi_aux.idx(row - 1, col), aux + fi_aux.idx(row, col), aux_global, 1,
         amdq + efi_asdq.down_edge(row, col), apdq + efi_asdq.down_edge(row, col),
         wave + efi_wave.down_edge(row, col), wave_speed + efi_s.down_edge(row, col));
    }

    for(row = r.begin(); row != std::min(r.end(), efi_wave.num_row_edge_transverse()); ++row){
      col = efi_wave.num_col_edge_transverse();
      rp(q + fi.idx(row, col - 1), q + fi.idx(row, col),
         aux + fi_aux.idx(row, col - 1), aux + fi_aux.idx(row, col), aux_global, 0,
         amdq + efi_asdq.left_edge(row, col), apdq + efi_asdq.left_edge(row, col),
         wave + efi_wave.left_edge(row, col), wave_speed + efi_s.left_edge(row, col));
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
  void_grid_eval_tbb_boundary boundary(rp, rp_grid_params, 
				       q, aux, aux_global, nx, ny, amdq, apdq, wave, wave_speed);

  tbb::parallel_for(::tbb::blocked_range2d<unsigned,unsigned>(1, efi.num_row_edge_transverse(),
							      1, efi.num_col_edge_transverse()),
		    body);

  tbb::parallel_for(::tbb::blocked_range<unsigned>(1, std::max(efi.num_row_edge_transverse(), efi.num_col_edge_transverse())), 
		    boundary);

}
