#include "benchmark_utils.h"

int main(int argc, char **argv)
{
  return benchmark(argc, argv,
		   num_euler_rp_grid_eval_kernels,
		   euler_rp_grid_eval_names,
		   euler_rp_grid_evals,
		   euler_rp_grid_params, 
		   &euler_rp_aux_global_default
		   );
}

