#include "euler_grid_eval.h"
#include "template_grid_eval.ipp"
#include "void_grid_eval.h"

const char * euler_rp_grid_eval_names[] =
  {
#ifdef  USE_TEMPLATE_GRID_EVAL
    "template",
#endif // USE_TEMPLATE_GRID_EVAL
    "void_star",
    "void_star_tbb",
    "void_star_omp",
  };

const rp_grid_eval_t euler_rp_grid_evals[] =
  {
#ifdef  USE_TEMPLATE_GRID_EVAL
    euler_rp_grid_eval_template,
#endif // USE_TEMPLATE_GRID_EVAL
    euler_rp_grid_eval_void_serial,
    euler_rp_grid_eval_void_tbb,
    euler_rp_grid_eval_void_omp,
    // TODO add other euler_rp_grid_eval functions here
  };

const size_t num_euler_rp_grid_eval_kernels = sizeof(euler_rp_grid_evals)/sizeof(rp_grid_eval_t);

#ifdef  USE_TEMPLATE_GRID_EVAL
void euler_rp_grid_eval_template(const real* q,  const real* aux, const void* aux_global,
				 const int nx, const  int ny,
				 real* amdq, real* apdq, real* wave,
				 real* wave_speed)
{
  template_rp_grid_eval_serial<euler_rp_aux_global_t>(&euler_rp,
                                                      euler_rp_grid_params,
                                                      q,
                                                      aux,
                                                      (euler_rp_aux_global_t*) aux_global,
                                                      nx,
                                                      ny,
                                                      amdq,
                                                      apdq,
                                                      wave,
                                                      wave_speed);
}
#endif // USE_TEMPLATE_GRID_EVAL

void euler_rp_grid_eval_void_serial(const real* q,  const real* aux, const void* aux_global,
				    const int nx, const  int ny,
				    real* amdq, real* apdq, real* wave,
				    real* wave_speed)
{
  void_rp_grid_eval_serial(&euler_rp,
                           euler_rp_grid_params,
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

void euler_rp_grid_eval_void_tbb(const real* q,  const real* aux,
                                     const void* aux_global,
                                     const int nx, const  int ny,
                                     real* amdq, real* apdq, real* wave,
                                     real* wave_speed)
{
  void_rp_grid_eval_tbb(&euler_rp,
                        euler_rp_grid_params,
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

void euler_rp_grid_eval_void_omp(const real* q,  const real* aux,
                                     const void* aux_global,
                                     const int nx, const  int ny,
                                     real* amdq, real* apdq, real* wave,
                                     real* wave_speed)
{
  void_rp_grid_eval_omp(&euler_rp,
                        euler_rp_grid_params,
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

