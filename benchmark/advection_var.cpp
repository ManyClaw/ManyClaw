#include "benchmark_utils.h"

int main(int argc, char **argv)
{
  return benchmark(argc, argv,
		   num_advection_var_rp_grid_eval_kernels,
		   advection_var_rp_grid_eval_names,
		   advection_var_rp_grid_evals,
		   advection_var_rp_grid_params, 
		   NULL
		   );
}
