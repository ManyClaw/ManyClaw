#include "euler_grid_eval.h"
#include "template_grid_eval.ipp"

const char * euler_rp_grid_eval_names[] =
  {
    "template"
  };

const rp_grid_eval_t euler_rp_grid_evals[] =
  {
    euler_rp_grid_eval_template
    // TODO add other euler_rp_grid_eval functions here
  };

const size_t num_euler_rp_grid_eval_kernels = sizeof(euler_rp_grid_evals)/sizeof(rp_grid_eval_t);


void euler_rp_grid_eval_template( const real* q,  const real* aux,
                                      const int nx, const  int ny,
                                      real* amdq, real* apdq, real* wave,
                                      real* wave_speed)
{
  template_rp_grid_eval_serial<euler_rp_aux_global_t>(&euler_rp,
                                                      euler_rp_grid_params,
                                                      euler_rp_aux_global,
                                                      q,
                                                      aux,
                                                      nx,
                                                      ny,
                                                      amdq,
                                                      apdq,
                                                      wave,
                                                      wave_speed);
}
