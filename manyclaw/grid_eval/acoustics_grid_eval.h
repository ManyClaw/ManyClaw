#ifndef ACOUSTICS_GRID_EVAL_H
#define ACOUSTICS_GRID_EVAL_H

#include "../common/data_structures.h"
#include "../riemann/acoustics_rp.h"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>

#include "template_grid_eval.h"

void acoustics_const_rp_grid_eval_template(const real* q,
					   const real* aux,
					   const void* aux_global,
					   const int nx,
					   const int ny,
					   real* amdq,
					   real* apdq,
					   real* wave,
					   real* wave_speed);

void acoustics_const_rp_grid_eval_void_serial(const real* q,  
					      const real* aux,
					      const void* aux_global,
					      const int nx, 
					      const int ny,
					      real* amdq, 
					      real* apdq, 
					      real* wave,
					      real* wave_speed);

void acoustics_const_rp_grid_eval_void_tbb(const real* q,  const real* aux, 
					   const void* aux_global,
					   const int nx, const  int ny,
					   real* amdq, real* apdq, real* wave,
					   real* wave_speed);

void acoustics_const_rp_grid_eval_void_omp(const real* q,  const real* aux, 
					   const void* aux_global,
					   const int nx, const  int ny,
					   real* amdq, real* apdq, real* wave,
					   real* wave_speed);

void acoustics_var_rp_grid_eval_template(const real* q,
					 const real* aux,
					 const void* aux_global,
					 const int nx,
					 const int ny,
					 real* amdq,
					 real* apdq,
					 real* wave,
					 real* wave_speed);

void acoustics_var_rp_grid_eval_void_serial(const real* q,  const real* aux, 
                                            const void* aux_global,
                                            const int nx, const  int ny,
                                            real* amdq, real* apdq, real* wave,
                                            real* wave_speed);

void acoustics_var_rp_grid_eval_void_tbb(const real* q,  const real* aux, 
                                         const void* aux_global,
                                         const int nx, const  int ny,
                                         real* amdq, real* apdq, real* wave,
                                         real* wave_speed);

void acoustics_var_rp_grid_eval_void_omp(const real* q,  const real* aux, 
                                         const void* aux_global,
                                         const int nx, const  int ny,
                                         real* amdq, real* apdq, real* wave,
                                         real* wave_speed);

extern const char * acoustics_const_rp_grid_eval_names[];
extern const rp_grid_eval_t acoustics_const_rp_grid_evals[];
extern const size_t num_acoustics_const_rp_grid_eval_kernels;

extern const char * acoustics_var_rp_grid_eval_names[];
extern const rp_grid_eval_t acoustics_var_rp_grid_evals[];
extern const size_t num_acoustics_var_rp_grid_eval_kernels;

#endif // ACOUSTICS_GRID_EVAL_H

