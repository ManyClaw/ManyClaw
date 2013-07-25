#include "advection_grid_eval.h"
#include "template_grid_eval.ipp"

const char * advection_rp_grid_eval_names[] =
  {
    "serial",
    "TBB",
    "omp",
    "template"
  };

const rp_grid_eval_t advection_rp_grid_evals[] =
  {
    advection_rp_grid_eval_serial,
    advection_rp_grid_eval_tbb,
    advection_rp_grid_eval_omp,
    advection_rp_grid_eval_template
    // TODO add other advection_rp_grid_eval functions here
  };

const size_t num_advection_rp_grid_eval_kernels = sizeof(advection_rp_grid_evals)/sizeof(rp_grid_eval_t);


void advection_rp_grid_eval_serial( const real* q,  const real* aux,
				    const int nx, const int ny,
				    real* amdq, real* apdq, real* wave,
				    real* wave_speed)
{
  int col, row;
  const int num_ghost = advection_rp_grid_params.num_ghost;
  const int num_eqn = advection_rp_grid_params.num_eqn;
  const int num_wave = advection_rp_grid_params.num_wave;

  FieldIndexer fi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_wave);

  for(row = 1; row < efi.num_row_edge_transverse(); ++row){
    for(col = 1; col < efi.num_col_edge_normal() + 1; ++col) {
      // printf("-> Calling %d, %d ", row, col);
      // printf("-> fi.idx %d, efi left %d down %d\n", fi.idx(row, col), efi.left_edge(row, col), efi.down_edge(row, col));
      advection_rp(q + fi.idx(row, col-1), q + fi.idx(row, col),
                   aux, aux,  &advection_rp_aux_global, 0,
                   amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
                   wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));

      advection_rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
                   aux, aux,  &advection_rp_aux_global, 1,
                   amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
		   wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
    }
  }

  for(col = 1; col < efi.num_col_edge_transverse(); ++col) {
    row = efi.num_row_edge_transverse();
    // printf("-> Calling %d, %d ", row, col);
    // printf("-> fi.idx %d, efi left %d down %d\n", fi.idx(row, col), efi.left_edge(row, col), efi.down_edge(row, col));
    advection_rp(q + fi.idx(row - 1, col), q + fi.idx(row, col),
                 aux, aux,  &advection_rp_aux_global, 1,
                 amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
                 wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
  }

  for(row = 1; row < efi.num_row_edge_transverse(); ++row){
    col = efi.num_col_edge_transverse();
    // printf("-> Calling %d, %d ", row, col);
    // printf("-> fi.idx %d, efi left %d down %d\n", fi.idx(row, col), efi.left_edge(row, col), efi.down_edge(row, col));
    advection_rp(q + fi.idx(row, col - 1), q + fi.idx(row, col),
                 aux, aux,  &advection_rp_aux_global, 0,
                 amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
                 wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));
  }
}


void advection_rp_grid_eval_omp( const real* q,  const real* aux,
                                 const int nx, const  int ny,
                                 real* amdq, real* apdq, real* wave,
                                 real* wave_speed)
{
  int col, row;
  const int num_ghost = advection_rp_grid_params.num_ghost;
  const int num_eqn = advection_rp_grid_params.num_eqn;
  const int num_wave = advection_rp_grid_params.num_wave;

  FieldIndexer fi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_wave);

#pragma omp parallel shared(q, aux, amdq, apdq, wave, wave_speed) private(col, row)
  {
#pragma omp for schedule(runtime) nowait
    for(row = 1; row < efi.num_row_edge_transverse(); ++row){
      for(col = 1; col < efi.num_col_edge_normal() + 1; ++col) {
	advection_rp(q + fi.idx(row, col-1), q + fi.idx(row, col),
		     aux, aux,  &advection_rp_aux_global, 0,
		     amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
		     wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));
	
        advection_rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
		     aux, aux,  &advection_rp_aux_global, 1,
		     amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
                     wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
      }
    }
    
#pragma omp for schedule(runtime) nowait
    for(col = 1; col < efi.num_col_edge_transverse(); ++col) {
      row = efi.num_row_edge_transverse();
      advection_rp(q + fi.idx(row - 1, col), q + fi.idx(row, col),
		   aux, aux,  &advection_rp_aux_global, 1,
		   amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
		   wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
    }
    
#pragma omp for schedule(runtime) nowait
    for(row = 1; row < efi.num_row_edge_transverse(); ++row){
      col = efi.num_col_edge_transverse();
      advection_rp(q + fi.idx(row, col - 1), q + fi.idx(row, col),
		   aux, aux,  &advection_rp_aux_global, 0,
		   amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
		   wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));
    }
  }
}


struct advection_rp_grid_eval_tbb_body
{
  const real* q;
  const real* aux;
  const int nx;
  const int ny;
  real* amdq; real* apdq;
  real* wave; real* wave_speed;

  advection_rp_grid_eval_tbb_body(const real* q,  const real* aux,
                             const int nx, const int ny,
                             real* amdq, real* apdq, real* wave,
                             real* wave_speed)
    : q(q), aux(aux), nx(nx), ny(ny), amdq(amdq), apdq(apdq), wave(wave), wave_speed(wave_speed)
  {}

  void operator()(const tbb::blocked_range2d<int>& r) const
  {

    int col, row;
    const int num_ghost = advection_rp_grid_params.num_ghost;
    const int num_eqn = advection_rp_grid_params.num_eqn;
    const int num_wave = advection_rp_grid_params.num_wave;

    FieldIndexer fi(nx, ny, num_ghost, num_eqn);
    EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_wave);

    for(row = r.rows().begin(); row < r.rows().end(); ++row){
      for(col = r.cols().begin(); col < r.cols().end() - 1; ++col) {
	advection_rp(q + fi.idx(row, col-1), q + fi.idx(row, col),
		     aux, aux,  &advection_rp_aux_global, 0,
		     amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
		     wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));
        advection_rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
		     aux, aux,  &advection_rp_aux_global, 1,
		     amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
                     wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
      }
    }

    for(col = r.cols().begin(); col < r.cols().end(); ++col) {
      row = r.rows().end();
      advection_rp(q + fi.idx(row - 1, col), q + fi.idx(row, col),
		   aux, aux,  &advection_rp_aux_global, 1,
		   amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
		   wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
    }

    for(row = r.rows().begin(); row < r.rows().end(); ++row){
      col = r.cols().end() + 1;
      advection_rp(q + fi.idx(row, col - 1), q + fi.idx(row, col),
		   aux, aux,  &advection_rp_aux_global, 0,
		   amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
		   wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));
    } 
  }
};


void advection_rp_grid_eval_tbb( const real* q,  const real* aux,
				 const int nx, const  int ny,
				 real* amdq, real* apdq, real* wave,
				 real* wave_speed)
{
  const int num_ghost = advection_rp_grid_params.num_ghost;
  const int num_eqn = advection_rp_grid_params.num_eqn;
  const int num_wave = advection_rp_grid_params.num_wave;

  FieldIndexer fi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_wave);
  advection_rp_grid_eval_tbb_body body(q, aux, nx, ny, amdq, apdq, wave, wave_speed);

  // note: we use nx+1 and ny+1 here and < in the body (instead of <= in the serial reference)
  tbb::parallel_for(::tbb::blocked_range2d<int,int>(1, efi.num_row_edge_transverse(), 
						    1, efi.num_col_edge_transverse()), body);
}

void advection_rp_grid_eval_template( const real* q,  const real* aux,
                                      const int nx, const  int ny,
                                      real* amdq, real* apdq, real* wave,
                                      real* wave_speed)
{
  template_rp_grid_eval_serial<advection_rp_aux_global_t>(&advection_rp,
                                                          advection_rp_grid_params,
                                                          advection_rp_aux_global,
                                                          q,
                                                          aux,
                                                          nx,
                                                          ny,
                                                          amdq,
                                                          apdq,
                                                          wave,
                                                          wave_speed);
}
