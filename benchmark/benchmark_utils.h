#ifndef __BENCHMARK_UTILS_H
#define __BENCHMARK_UTILS_H

#include <manyclaw/manyclaw.h>
#include <tbb/task_scheduler_init.h>

#include "timer.h"

double benchmark_grid_eval(int nx, int ny, rp_grid_params_t params,
                           rp_grid_eval_t rp_grid_eval,
                           const void* aux_global);

int benchmark(int argc, char **argv, size_t num_kernel, const char *names[], const rp_grid_eval_t grid_evals[],
	      rp_grid_params_t params, const void* aux_global);

#endif // __BENCHMARK_UTILS_H
