#ifndef ADVECTION_GRID_EVAL_H
#define ADVECTION_GRID_EVAL_H

#include "../common/data_structures.h"
#include "../riemann/advection_rp.h"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>

#include "template_grid_eval.h"

// Evaluates advection_rp via serial execution
void advection_rp_grid_eval_serial(const real* q,
                   const real* aux,
                                   const void* aux_global,
                   const int nx,
                   const int ny,
                   real* amdq,
                   real* apdq,
                   real* wave,
                   real* wave_speed);

// Evaluates advection_rp via parallel execution via openmp
void advection_rp_grid_eval_omp(const real* q,
                const real* aux,
                                const void* aux_global,
                const int nx,
                const int ny,
                real* amdq,
                real* apdq,
                real* wave,
                real* wave_speed);

// Evaluates advection_rp via parallel execution via TBB
void advection_rp_grid_eval_tbb(const real* q,
                const real* aux,
                                const void* aux_global,
                const int nx,
                const int ny,
                real* amdq,
                real* apdq,
                real* wave,
                real* wave_speed);

void advection_rp_grid_eval_template(const real* q,
                                     const real* aux,
                                     const void* aux_global,
                                     const int nx,
                                     const int ny,
                                     real* amdq,
                                     real* apdq,
                                     real* wave,
                                     real* wave_speed);

void advection_rp_grid_eval_void( const real* q,  
                                  const real* aux,
                                  const void* aux_global,
                                  const int nx, 
                                  const int ny,
                                  real* amdq, 
                                  real* apdq, 
                                  real* wave,
                                  real* wave_speed);

void advection_var_rp_grid_eval_void( const real* q,  const real* aux, 
                                  const void* aux_global,
                                  const int nx, const  int ny,
                                  real* amdq, real* apdq, real* wave,
                                  real* wave_speed);

extern const char * advection_rp_grid_eval_names[];
extern const rp_grid_eval_t advection_rp_grid_evals[];
extern const size_t num_advection_rp_grid_eval_kernels;

#endif // ADVECTION_GRID_EVAL_H

