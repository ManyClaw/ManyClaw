#include <limits.h>
#include <manyclaw/manyclaw.h>
#include "gtest/gtest.h"

::testing::AssertionResult ArraysMatch(const real *expected,
                                       const real *actual, int size){
  std::vector<int> bad_loc;
  for (int i=0; i < size; ++i){
    if (std::fabs(expected[i] - actual[i]) > 1e-8){
      bad_loc.push_back(i);
    }
  }

  if (bad_loc.size() > 0) {
    std::ostringstream fail;
    fail << std::endl;
    for(size_t i=0; i < bad_loc.size(); ++i){
      int idx = bad_loc[i];
      fail << std::scientific
           << "  array[" << idx
           << "] (" << actual[idx] << ") != expected[" << idx
           << "] (" << expected[idx] << ")\n";
    }
    fail << "  Num_wrong:" << bad_loc.size();
    return ::testing::AssertionFailure() << fail.str();
  }

  return ::testing::AssertionSuccess();
}

// Test the indexer methods
TEST(AdvectionStepper, base) {
  int nx = 5, ny = 5;
  int num_ghost = advection_rp_grid_params.num_ghost; 
  int num_eqns = advection_rp_grid_params.num_eqn; 
  int num_waves = advection_rp_grid_params.num_waves; 
  FieldIndexer fi(nx, ny, num_ghost, num_eqns);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqns, num_waves);
  real q[fi.size()];
  real amdq[efi.size()], amdq_gold[efi.size()];
  real apdq[efi.size()], apdq_gold[efi.size()];
  real wave[efi.size()], wave_gold[efi.size()];
  real speed[efi.size()], speed_gold[efi.size()];

  for(int idx=0; idx < fi.row_size() * fi.col_size(); ++idx)
      q[idx] = 0.0;
  q[fi.idx(5, 5)] = 1.0;

  for(int idx=0; idx < fi.row_size() * fi.col_size(); ++idx) {
    wave_gold[idx] = 0.0;
    speed_gold[idx] = 1.0;
    amdq_gold[idx] = 0.0;
    apdq_gold[idx] = 0.0;
  }

  wave_gold[efi.left_edge(5,5)] = 1.0;
  speed_gold[efi.left_edge(5,5)] = 1.0;
  amdq_gold[efi.left_edge(5,5)] = 0.0;
  apdq_gold[efi.left_edge(5,5)] = 1.0;

  wave_gold[efi.right_edge(5,5)] = -1.0;
  speed_gold[efi.right_edge(5,5)] = 1.0;
  amdq_gold[efi.right_edge(5,5)] = 0.0;
  apdq_gold[efi.right_edge(5,5)] = -1.0;

  wave_gold[efi.down_edge(5,5)] = 1.0;
  speed_gold[efi.down_edge(5,5)] = 1.0;
  amdq_gold[efi.down_edge(5,5)] = 0.0;
  apdq_gold[efi.down_edge(5,5)] = 1.0;

  wave_gold[efi.up_edge(5,5)] = -1.0;
  speed_gold[efi.up_edge(5,5)] = 1.0;
  amdq_gold[efi.up_edge(5,5)] = 0.0;
  apdq_gold[efi.up_edge(5,5)] = -1.0;
  
  advection_rp_step_serial_cellwise(q, NULL, nx, ny, amdq, apdq, wave, speed);

  EXPECT_TRUE(ArraysMatch(wave, wave_gold, efi.size()));
  EXPECT_TRUE(ArraysMatch(speed, speed_gold, efi.size()));
  EXPECT_TRUE(ArraysMatch(amdq, amdq_gold, efi.size()));
  EXPECT_TRUE(ArraysMatch(apdq, apdq_gold, efi.size()));
}
