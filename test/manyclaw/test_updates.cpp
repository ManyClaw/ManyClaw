#include "gtest/gtest.h"
#include "../test_utils.h"

#include <manyclaw/manyclaw.h>

#include <limits.h>

// Test the indexer methods
TEST(UpdateTest, FirstOrderDimensionalSplitting) {
  int nx = 3, ny = 3;
  int num_ghost = advection_rp_grid_params.num_ghost; 
  int num_eqn = advection_rp_grid_params.num_eqn; 
  int num_wave = advection_rp_grid_params.num_wave; 
  real dtdx = 1.0;
  FieldIndexer fi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn, num_wave);
  real q[fi.size()], q_gold[fi.size()];
  real amdq[efi.size()];
  real apdq[efi.size()];
  real wave[efi.size()];
  real speed[efi.size()];

  for(unsigned idx=0; idx < efi.size(); ++idx){
    amdq[idx] = 0.0;
    apdq[idx] = 0.0;
    wave[idx] = 0.0;
    speed[idx] = 0.0;
  }

  int r=1, c=1;
  amdq[efi.down_edge(r, c)] = -1.0;
  apdq[efi.left_edge(r, c)] = 1.0;

  for(unsigned idx=0; idx < fi.size(); ++idx) {
    q[idx] = 0.0;
    q_gold[idx] = 0.0;
  }
  q_gold[fi.idx(r, c)] = -1.0;
  q_gold[fi.down(r, c)] = 1.0;

  updater_first_order_dimensional_splitting(q, NULL, nx, ny, amdq, apdq, wave, speed, num_ghost, num_eqn, dtdx);

  EXPECT_TRUE(ArraysMatch(q, q_gold, fi.size()));
}
