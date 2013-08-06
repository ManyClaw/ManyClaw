#include "benchmark_utils.h"

int main(int argc, char **argv)
{
  return benchmark(argc, argv,
		   num_acoustics_var_rp_grid_eval_kernels,
		   acoustics_var_rp_grid_eval_names,
		   acoustics_var_rp_grid_evals,
		   acoustics_var_rp_grid_params, 
		   &acoustics_var_rp_aux_global_default
		   );
}

