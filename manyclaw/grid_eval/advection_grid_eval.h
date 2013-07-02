#ifndef ADVECTION_GRID_EVAL_H
#define ADVECTION_GRID_EVAL_H

#include "../common/data_structures.h"
#include "../riemann/advection_rp.h"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>


// Evaluates advection_rp via serial execution
void advection_rp_grid_eval_serial(const real* q,
				   const real* aux,
				   const int nx,
				   const int ny,
				   real* amdq,
				   real* apdq,
				   real* wave,
				   real* wave_speed);

// Evaluates advection_rp via parallel execution via openmp
void advection_rp_grid_eval_omp(const real* q,
				const real* aux,
				const int nx,
				const int ny,
				real* amdq,
				real* apdq,
				real* wave,
				real* wave_speed);

// Evaluates advection_rp via parallel execution via TBB
void advection_rp_grid_eval_tbb(const real* q,
				const real* aux,
				const int nx,
				const int ny,
				real* amdq,
				real* apdq,
				real* wave,
				real* wave_speed);


const char * advection_rp_grid_eval_names[] =
  {
    "serial",
    "TBB",
    "omp"
  };

rp_grid_eval_t advection_rp_grid_evals[] =
  {
    advection_rp_grid_eval_serial,
    advection_rp_grid_eval_tbb,
    advection_rp_grid_eval_omp
    // TODO add other advection_rp_grid_eval functions here
  };

size_t num_advection_rp_grid_eval_kernels = sizeof(advection_rp_grid_evals)/sizeof(rp_grid_eval_t);

#endif // ADVECTION_GRID_EVAL_H

