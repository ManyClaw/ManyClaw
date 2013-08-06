#include "gtest/gtest.h"
#include "../../test_utils.h"

#include <manyclaw/manyclaw.h>

#include <limits.h>

// validates the different methods
TEST(EulerGridEval, DISABLED_validate) {
  int nx = 10, ny = 10;
  rp_grid_params_t params = euler_rp_grid_params; 
  const euler_rp_aux_global_t aux_global = euler_rp_aux_global_default;

  for (unsigned idx = 1; idx < num_euler_rp_grid_eval_kernels; ++idx){
    EXPECT_TRUE(GridEvalsMatch(nx, ny, params, 
			       euler_rp_grid_evals[idx], euler_rp_grid_eval_names[idx], 
			       euler_rp_grid_evals[0], euler_rp_grid_eval_names[0], 
			       &aux_global) );
  }
  

}

