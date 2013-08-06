#include "gtest/gtest.h"
#include "../../test_utils.h"

#include <manyclaw/manyclaw.h>

#include <limits.h>

// validates the different methods
TEST(AcousticsGridEval, DISABLED_validate_const) {
  int nx = 10, ny = 10;
  rp_grid_params_t params = acoustics_const_rp_grid_params; 
  const acoustics_const_rp_aux_global_t aux_global =  acoustics_const_rp_aux_global_default;

  for (unsigned idx = 1; idx < num_acoustics_const_rp_grid_eval_kernels; ++idx){
    EXPECT_TRUE(GridEvalsMatch(nx, ny, params, 
			       acoustics_const_rp_grid_evals[idx], acoustics_const_rp_grid_eval_names[idx], 
			       acoustics_const_rp_grid_evals[0], acoustics_const_rp_grid_eval_names[0], 
			       &aux_global) );
  }
  

}

TEST(AcousticsGridEval, DISABLED_validate_var) {
  int nx = 10, ny = 10;
  rp_grid_params_t params = acoustics_var_rp_grid_params; 
  const acoustics_var_rp_aux_global_t aux_global =  acoustics_var_rp_aux_global_default;

  for (unsigned idx = 1; idx < num_acoustics_var_rp_grid_eval_kernels; ++idx){
    EXPECT_TRUE(GridEvalsMatch(nx, ny, params, 
			       acoustics_var_rp_grid_evals[idx], acoustics_var_rp_grid_eval_names[idx], 
			       acoustics_var_rp_grid_evals[0], acoustics_var_rp_grid_eval_names[0], 
			       &aux_global) );
  }
  

}
