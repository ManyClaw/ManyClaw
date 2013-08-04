#include "advection_grid_eval.h"
#include "template_grid_eval.ipp"
#include "void_grid_eval.h"

const char * advection_rp_grid_eval_names[] =
  {
    "serial",
    "TBB",
    "omp",
#ifdef  USE_TEMPLATE_GRID_EVAL
    "template",
#endif // USE_TEMPLATE_GRID_EVAL
    "void_star",
    "void_star_tbb",
    "void_star_omp",
  };

const rp_grid_eval_t advection_rp_grid_evals[] =
  {
    advection_rp_grid_eval_serial,
    advection_rp_grid_eval_tbb,
    advection_rp_grid_eval_omp,
#ifdef  USE_TEMPLATE_GRID_EVAL
    advection_rp_grid_eval_template,
#endif // USE_TEMPLATE_GRID_EVAL
    advection_rp_grid_eval_void_serial,
    advection_rp_grid_eval_void_tbb,
    advection_rp_grid_eval_void_omp,
    // TODO add other advection_rp_grid_eval functions here
  };


const size_t num_advection_rp_grid_eval_kernels = sizeof(advection_rp_grid_evals)/sizeof(rp_grid_eval_t);

const char * advection_var_rp_grid_eval_names[] =
  {
    "void_star",
    "void_star_tbb",
    "void_star_omp",
  };

const rp_grid_eval_t advection_var_rp_grid_evals[] =
  {
    advection_var_rp_grid_eval_void_serial,
    advection_var_rp_grid_eval_void_tbb,
    advection_var_rp_grid_eval_void_omp,
    // TODO add other advection_rp_grid_eval functions here
  };

const size_t num_advection_var_rp_grid_eval_kernels = sizeof(advection_var_rp_grid_evals)/sizeof(rp_grid_eval_t);


void advection_rp_grid_eval_serial( const real* q,  const real* aux, const void* aux_global,
                                    const int nx, const int ny,
                                    real* amdq, real* apdq, real* wave,
                                    real* wave_speed)
{
  unsigned col, row;
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
                   aux, aux,  aux_global, 0,
                   amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
                   wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));

      advection_rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
                   aux, aux,  aux_global, 1,
                   amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
                   wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
    }
  }

  for(col = 1; col < efi.num_col_edge_transverse(); ++col) {
    row = efi.num_row_edge_transverse();
    // printf("-> Calling %d, %d ", row, col);
    // printf("-> fi.idx %d, efi left %d down %d\n", fi.idx(row, col), efi.left_edge(row, col), efi.down_edge(row, col));
    advection_rp(q + fi.idx(row - 1, col), q + fi.idx(row, col),
                 aux, aux,  aux_global, 1,
                 amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
                 wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
  }

  for(row = 1; row < efi.num_row_edge_transverse(); ++row){
    col = efi.num_col_edge_transverse();
    // printf("-> Calling %d, %d ", row, col);
    // printf("-> fi.idx %d, efi left %d down %d\n", fi.idx(row, col), efi.left_edge(row, col), efi.down_edge(row, col));
    advection_rp(q + fi.idx(row, col - 1), q + fi.idx(row, col),
                 aux, aux,  aux_global, 0,
                 amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
                 wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));
  }
}


void advection_rp_grid_eval_omp( const real* q,  const real* aux, const void* aux_global,
                                 const int nx, const  int ny,
                                 real* amdq, real* apdq, real* wave,
                                 real* wave_speed)
{
  unsigned col, row;
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
                     aux, aux,  aux_global, 0,
                     amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
                     wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));

        advection_rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
                     aux, aux,  aux_global, 1,
                     amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
                     wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
      }
    }

#pragma omp for schedule(runtime) nowait
    for(col = 1; col < efi.num_col_edge_transverse(); ++col) {
      row = efi.num_row_edge_transverse();
      advection_rp(q + fi.idx(row - 1, col), q + fi.idx(row, col),
                   aux, aux,  aux_global, 1,
                   amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
                   wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
    }

#pragma omp for schedule(runtime) nowait
    for(row = 1; row < efi.num_row_edge_transverse(); ++row){
      col = efi.num_col_edge_transverse();
      advection_rp(q + fi.idx(row, col - 1), q + fi.idx(row, col),
                   aux, aux,  aux_global, 0,
                   amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
                   wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));
    }
  }
}


struct advection_rp_grid_eval_tbb_body
{
  const real* q;
  const real* aux;
  const void* aux_global;
  const int nx;
  const int ny;
  real* amdq; real* apdq;
  real* wave; real* wave_speed;

  advection_rp_grid_eval_tbb_body(const real* q,  const real* aux, const void* aux_global,
                             const int nx, const int ny,
                             real* amdq, real* apdq, real* wave,
                             real* wave_speed)
    : q(q), aux(aux), aux_global(aux_global), nx(nx), ny(ny), amdq(amdq), apdq(apdq), wave(wave), wave_speed(wave_speed)
  {}

  void operator()(const tbb::blocked_range2d<unsigned>& r) const
  {

    unsigned col, row;
    const int num_ghost = advection_rp_grid_params.num_ghost;
    const int num_eqn = advection_rp_grid_params.num_eqn;
    const int num_wave = advection_rp_grid_params.num_wave;

    FieldIndexer fi(nx, ny, num_ghost, num_eqn);
    EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_wave);

    for(row = r.rows().begin(); row != r.rows().end(); ++row){
      for(col = r.cols().begin(); col != std::min(r.cols().end(), efi.num_col_edge_normal() + 1); ++col) {
        advection_rp(q + fi.idx(row, col-1), q + fi.idx(row, col),
                     aux, aux,  aux_global, 0,
                     amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
                     wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));
        advection_rp(q + fi.idx(row-1, col), q + fi.idx(row, col),
                     aux, aux,  aux_global, 1,
                     amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
                     wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
      }
    }
  }
};

struct advection_rp_grid_eval_tbb_boundary
{
  const real* q;
  const real* aux;
  const void* aux_global;
  const int nx;
  const int ny;
  real* amdq; real* apdq;
  real* wave; real* wave_speed;

  advection_rp_grid_eval_tbb_boundary(const real* q,  const real* aux, const void* aux_global,
				      const int nx, const int ny,
				      real* amdq, real* apdq, real* wave,
				      real* wave_speed)
    : q(q), aux(aux), aux_global(aux_global), nx(nx), ny(ny), amdq(amdq), apdq(apdq), wave(wave), wave_speed(wave_speed)
  {}

  void operator()(const tbb::blocked_range<unsigned>& r) const
  {
    
    unsigned col, row;
    const int num_ghost = advection_rp_grid_params.num_ghost;
    const int num_eqn = advection_rp_grid_params.num_eqn;
    const int num_wave = advection_rp_grid_params.num_wave;

    FieldIndexer fi(nx, ny, num_ghost, num_eqn);
    EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_wave);
    
    for(col = r.begin(); col != std::min(r.end(), efi.num_col_edge_transverse()); ++col) {
      row = efi.num_row_edge_transverse();
      advection_rp(q + fi.idx(row - 1, col), q + fi.idx(row, col),
                   aux, aux,  aux_global, 1,
                   amdq + efi.down_edge(row, col), apdq + efi.down_edge(row, col),
                   wave + efi.down_edge(row, col), wave_speed + efi.down_edge(row, col));
    }
    
    for(row = r.begin(); row != std::min(r.end(), efi.num_row_edge_transverse()); ++row){
      col = efi.num_col_edge_transverse();
      advection_rp(q + fi.idx(row, col - 1), q + fi.idx(row, col),
                   aux, aux,  aux_global, 0,
                   amdq + efi.left_edge(row, col), apdq + efi.left_edge(row, col),
                   wave + efi.left_edge(row, col), wave_speed + efi.left_edge(row, col));
    }
  }	  
};


void advection_rp_grid_eval_tbb( const real* q,  const real* aux, const void* aux_global,
                                 const int nx, const  int ny,
                                 real* amdq, real* apdq, real* wave,
                                 real* wave_speed)
{
  const int num_ghost = advection_rp_grid_params.num_ghost;
  const int num_eqn = advection_rp_grid_params.num_eqn;
  const int num_wave = advection_rp_grid_params.num_wave;

  FieldIndexer fi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_wave);
  advection_rp_grid_eval_tbb_body body(q, aux, aux_global, nx, ny, amdq, apdq, wave, wave_speed);
  advection_rp_grid_eval_tbb_boundary boundary(q, aux, aux_global, nx, ny, amdq, apdq, wave, wave_speed);

  tbb::parallel_for(::tbb::blocked_range2d<unsigned, unsigned>(1, efi.num_row_edge_transverse(),
							       1, efi.num_col_edge_transverse()),
		    body);
  tbb::parallel_for(::tbb::blocked_range<unsigned>(1, std::max(efi.num_row_edge_transverse(), efi.num_col_edge_transverse())), 
		    boundary);
}

// Macro this off since ForestClaw has trouble with templates.
#ifdef USE_TEMPLATE_GRID_EVAL
void advection_rp_grid_eval_template( const real* q,  const real* aux, const void* aux_global,
                                      const int nx, const  int ny,
                                      real* amdq, real* apdq, real* wave,
                                      real* wave_speed)
{
  template_rp_grid_eval_serial<advection_rp_aux_global_t>(&advection_rp,
                                                          advection_rp_grid_params,
                                                          q,
                                                          aux,
                                                          (advection_rp_aux_global*) aux_global,
                                                          nx,
                                                          ny,
                                                          amdq,
                                                          apdq,
                                                          wave,
                                                          wave_speed);
}
#endif // USE_TEMPLATE_GRID_EVAL

void advection_rp_grid_eval_void_serial( const real* q,  const real* aux, const void* aux_global,
                                  const int nx, const  int ny,
                                  real* amdq, real* apdq, real* wave,
                                  real* wave_speed)
{
  void_rp_grid_eval_serial(&advection_rp,
                           advection_rp_grid_params,
                           q,
                           aux,
                           aux_global,
                           nx,
                           ny,
                           amdq,
                           apdq,
                           wave,
                           wave_speed);
}


void advection_rp_grid_eval_void_tbb(const real* q,  const real* aux,
                                     const void* aux_global,
                                     const int nx, const  int ny,
                                     real* amdq, real* apdq, real* wave,
                                     real* wave_speed)
{
  void_rp_grid_eval_tbb(&advection_rp,
                        advection_rp_grid_params,
                        q,
                        aux,
                        aux_global,
                        nx,
                        ny,
                        amdq,
                        apdq,
                        wave,
                        wave_speed);
}

void advection_rp_grid_eval_void_omp(const real* q,  const real* aux,
                                     const void* aux_global,
                                     const int nx, const  int ny,
                                     real* amdq, real* apdq, real* wave,
                                     real* wave_speed)
{
  void_rp_grid_eval_omp(&advection_rp,
                        advection_rp_grid_params,
                        q,
                        aux,
                        aux_global,
                        nx,
                        ny,  
                        amdq,
                        apdq,
                        wave,
                        wave_speed);
}



void advection_var_rp_grid_eval_void_serial(const real* q,  const real* aux, const void* aux_global,
                                            const int nx, const  int ny,
                                            real* amdq, real* apdq, real* wave,
                                            real* wave_speed)
{
  void_rp_grid_eval_serial(&advection_var_rp,
                           advection_var_rp_grid_params,
                           q,
                           aux,
                           NULL,
                           nx,
                           ny,
                           amdq,
                           apdq,
                           wave,
                           wave_speed);
}


void advection_var_rp_grid_eval_void_tbb(const real* q,  const real* aux,
                                         const void* aux_global,
                                         const int nx, const  int ny,
                                         real* amdq, real* apdq, real* wave,
                                         real* wave_speed)
{
  void_rp_grid_eval_tbb(&advection_var_rp,
                        advection_var_rp_grid_params,
                        q,
                        aux,
                        NULL,
                        nx,
                        ny,
                        amdq,
                        apdq,
                        wave,
                        wave_speed);
}

void advection_var_rp_grid_eval_void_omp(const real* q,  const real* aux,
                                         const void* aux_global,
                                         const int nx, const  int ny,
                                         real* amdq, real* apdq, real* wave,
                                         real* wave_speed)
{
  void_rp_grid_eval_omp(&advection_var_rp,
                        advection_var_rp_grid_params,
                        q,
                        aux,
                        NULL,
                        nx,
                        ny,  
                        amdq,
                        apdq,
                        wave,
                        wave_speed);
}

