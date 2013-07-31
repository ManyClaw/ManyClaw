#ifndef VOID_GRID_EVAL_H
#define VOID_GRID_EVAL_H

#include "../common/data_structures.h"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>


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
                              real* wave_speed);

void void_rp_grid_eval_tbb(rp_t rp,
                           rp_grid_params_t rp_grid_params,
                           const real* q,
                           const real* aux,
                           const void* aux_global,
                           const int nx,
                           const int ny,
                           real* amdq,
                           real* apdq,
                           real* wave,
                           real* wave_speed);

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
                           real* wave_speed);

#endif // VOID_GRID_EVAL_H
