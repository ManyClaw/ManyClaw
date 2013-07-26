#ifndef TEMPLATE_GRID_EVAL_H
#define TEMPLATE_GRID_EVAL_H

// Macro this off since ForestClaw has trouble with templates.
#ifdef USE_TEMPLATE_GRID_EVAL

#include "../common/data_structures.h"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>


// Evaluates template_rp via serial execution
template<class rp_aux_global_t>
void template_rp_grid_eval_serial(rp_t rp,
                                  rp_grid_params_t rp_grid_params,
                                  rp_aux_global_t rp_aux_global,
                                  const real* q,
                                  const real* aux,
                                  const int nx,
                                  const int ny,
                                  real* amdq,
                                  real* apdq,
                                  real* wave,
                                  real* wave_speed);
#endif // USE_TEMPLATE_GRID_EVAL
#endif // TEMPLATE_GRID_EVAL_H
