#ifndef EULER_GRID_EVAL_H
#define EULER_GRID_EVAL_H

#include "../common/data_structures.h"
#include "../riemann/euler_rp.h"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>

#include "template_grid_eval.h"

void euler_rp_grid_eval_template(const real* q,
                                     const real* aux,
                                     const int nx,
                                     const int ny,
                                     real* amdq,
                                     real* apdq,
                                     real* wave,
                                     real* wave_speed);

extern const char * euler_rp_grid_eval_names[];
extern const rp_grid_eval_t euler_rp_grid_evals[];
extern const size_t num_euler_rp_grid_eval_kernels;

#endif // EULER_GRID_EVAL_H

