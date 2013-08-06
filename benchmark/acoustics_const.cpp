#include "benchmark_utils.h"

int main(int argc, char **argv)
{
  return benchmark(argc, argv,
		   num_acoustics_const_rp_grid_eval_kernels,
		   acoustics_const_rp_grid_eval_names,
		   acoustics_const_rp_grid_evals,
		   acoustics_const_rp_grid_params, 
		   &acoustics_const_rp_aux_global_default
		   );
}

