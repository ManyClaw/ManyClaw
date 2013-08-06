#include "acoustics_grid_eval.h"
#include "template_grid_eval.ipp"
#include "void_grid_eval.h"

const char * acoustics_const_rp_grid_eval_names[] =
  {
#ifdef  USE_TEMPLATE_GRID_EVAL
    "template",
#endif // USE_TEMPLATE_GRID_EVAL
    "void_star",
    "void_star_tbb",
    "void_star_omp",
  };

const rp_grid_eval_t acoustics_const_rp_grid_evals[] =
  {
#ifdef  USE_TEMPLATE_GRID_EVAL
    acoustics_const_rp_grid_eval_template,
#endif // USE_TEMPLATE_GRID_EVAL
    acoustics_const_rp_grid_eval_void_serial,
    acoustics_const_rp_grid_eval_void_tbb,
    acoustics_const_rp_grid_eval_void_omp,
   // TODO add other acoustics_const_rp_grid_eval functions here
  };

const size_t num_acoustics_const_rp_grid_eval_kernels = sizeof(acoustics_const_rp_grid_evals)/sizeof(rp_grid_eval_t);

const char * acoustics_var_rp_grid_eval_names[] =
  {
#ifdef  USE_TEMPLATE_GRID_EVAL
    "template",
#endif // USE_TEMPLATE_GRID_EVAL
    "void_star",
    "void_star_tbb",
    "void_star_omp",
  };

const rp_grid_eval_t acoustics_var_rp_grid_evals[] =
  {
#ifdef  USE_TEMPLATE_GRID_EVAL
    acoustics_var_rp_grid_eval_template,
#endif // USE_TEMPLATE_GRID_EVAL
    acoustics_var_rp_grid_eval_void_serial,
    acoustics_var_rp_grid_eval_void_tbb,
    acoustics_var_rp_grid_eval_void_omp,
   // TODO add other acoustics_var_rp_grid_eval functions here
  };

const size_t num_acoustics_var_rp_grid_eval_kernels = sizeof(acoustics_var_rp_grid_evals)/sizeof(rp_grid_eval_t);

#ifdef  USE_TEMPLATE_GRID_EVAL
void acoustics_const_rp_grid_eval_template( const real* q,  const real* aux, const void* aux_global,
					    const int nx, const  int ny,
					    real* amdq, real* apdq, real* wave,
					    real* wave_speed)
{
  template_rp_grid_eval_serial<acoustics_const_rp_aux_global_t>(&acoustics_const_rp,
                                                      acoustics_const_rp_grid_params,
                                                      q,
                                                      aux,
                                                      (acoustics_const_rp_aux_global_t*) aux_global,
                                                      nx,
                                                      ny,
                                                      amdq,
                                                      apdq,
                                                      wave,
                                                      wave_speed);
}
#endif // USE_TEMPLATE_GRID_EVAL

void acoustics_const_rp_grid_eval_void_serial(const real* q,  const real* aux, const void* aux_global,
					      const int nx, const  int ny,
					      real* amdq, real* apdq, real* wave,
					      real* wave_speed)
{
  void_rp_grid_eval_serial(&acoustics_const_rp,
                           acoustics_const_rp_grid_params,
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


void acoustics_const_rp_grid_eval_void_tbb(const real* q,  const real* aux,
					   const void* aux_global,
					   const int nx, const  int ny,
					   real* amdq, real* apdq, real* wave,
					   real* wave_speed)
{
  void_rp_grid_eval_tbb(&acoustics_const_rp,
                        acoustics_const_rp_grid_params,
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

void acoustics_const_rp_grid_eval_void_omp(const real* q,  const real* aux,
					   const void* aux_global,
					   const int nx, const  int ny,
					   real* amdq, real* apdq, real* wave,
					   real* wave_speed)
{
  void_rp_grid_eval_omp(&acoustics_const_rp,
                        acoustics_const_rp_grid_params,
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


#ifdef  USE_TEMPLATE_GRID_EVAL
void acoustics_var_rp_grid_eval_template( const real* q,  const real* aux, const void* aux_global,
                                      const int nx, const  int ny,
                                      real* amdq, real* apdq, real* wave,
                                      real* wave_speed)
{
  template_rp_grid_eval_serial<acoustics_var_rp_aux_global_t>(&acoustics_var_rp,
                                                      acoustics_var_rp_grid_params,
                                                      q,
                                                      aux,
                                                      (acoustics_var_rp_aux_global_t*) aux_global,
                                                      nx,
                                                      ny,
                                                      amdq,
                                                      apdq,
                                                      wave,
                                                      wave_speed);
}
#endif // USE_TEMPLATE_GRID_EVAL

void acoustics_var_rp_grid_eval_void_serial( const real* q,  const real* aux, const void* aux_global,
                              const int nx, const  int ny,
                              real* amdq, real* apdq, real* wave,
                              real* wave_speed)
{
  void_rp_grid_eval_serial(&acoustics_var_rp,
                           acoustics_var_rp_grid_params,
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


void acoustics_var_rp_grid_eval_void_tbb(const real* q,  const real* aux,
					 const void* aux_global,
					 const int nx, const  int ny,
					 real* amdq, real* apdq, real* wave,
                                     real* wave_speed)
{
  void_rp_grid_eval_tbb(&acoustics_var_rp,
                        acoustics_var_rp_grid_params,
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

void acoustics_var_rp_grid_eval_void_omp(const real* q,  const real* aux,
                                     const void* aux_global,
                                     const int nx, const  int ny,
                                     real* amdq, real* apdq, real* wave,
                                     real* wave_speed)
{
  void_rp_grid_eval_omp(&acoustics_var_rp,
                        acoustics_var_rp_grid_params,
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


