#include "benchmark_utils.h"

int main(int argc, char **argv)
{
  return benchmark(argc, argv,
		   num_advection_rp_grid_eval_kernels,
		   advection_rp_grid_eval_names,
		   advection_rp_grid_evals,
		   advection_rp_grid_params, 
		   &advection_rp_aux_global_default
		   );
}
