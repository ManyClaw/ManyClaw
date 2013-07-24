#include "gtest/gtest.h"
#include "../../test_utils.h"

#include <manyclaw/manyclaw.h>

#include <limits.h>

// Test the indexer methods
TEST(AdvectionStepper, base) {
  int nx = 3, ny = 3;
  int num_ghost = advection_rp_grid_params.num_ghost; 
  int num_eqns = advection_rp_grid_params.num_eqn; 
  int num_wave = advection_rp_grid_params.num_wave; 
  FieldIndexer fi(nx, ny, num_ghost, num_eqns);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqns, num_wave);
  real q[fi.size()];
  real amdq[efi.size()], amdq_gold[efi.size()];
  real apdq[efi.size()], apdq_gold[efi.size()];
  real wave[efi.size()], wave_gold[efi.size()];
  real speed[efi.size()], speed_gold[efi.size()];

  advection_rp_aux_global.u[0] = 1.0;
  advection_rp_aux_global.u[1] = 1.0;

  for(int idx=0; idx < fi.row_size() * fi.col_size(); ++idx)
      q[idx] = 0.0;
  int r=2, c=2;
  q[fi.idx(r, c)] = 1.0;

  for(int idx=0; idx < efi.size(); ++idx) {
    wave_gold[idx] = 0.0;
    speed_gold[idx] = 1.0;
    amdq_gold[idx] = 0.0;
    apdq_gold[idx] = 0.0;
  }

  wave_gold[efi.left_edge(r, c)] = 1.0;
  speed_gold[efi.left_edge(r, c)] = 1.0;
  amdq_gold[efi.left_edge(r, c)] = 0.0;
  apdq_gold[efi.left_edge(r, c)] = 1.0;

  wave_gold[efi.right_edge(r, c)] = -1.0;
  speed_gold[efi.right_edge(r, c)] = 1.0;
  amdq_gold[efi.right_edge(r, c)] = 0.0;
  apdq_gold[efi.right_edge(r, c)] = -1.0;

  wave_gold[efi.down_edge(r, c)] = 1.0;
  speed_gold[efi.down_edge(r, c)] = 1.0;
  amdq_gold[efi.down_edge(r, c)] = 0.0;
  apdq_gold[efi.down_edge(r, c)] = 1.0;

  wave_gold[efi.up_edge(r, c)] = -1.0;
  speed_gold[efi.up_edge(r, c)] = 1.0;
  amdq_gold[efi.up_edge(r, c)] = 0.0;
  apdq_gold[efi.up_edge(r, c)] = -1.0;
  
  advection_rp_grid_eval_serial(q, NULL, nx, ny, amdq, apdq, wave, speed);

  EXPECT_TRUE(ArraysMatch(wave, wave_gold, efi.size()));
  EXPECT_TRUE(ArraysMatch(speed, speed_gold, efi.size()));
  EXPECT_TRUE(ArraysMatch(amdq, amdq_gold, efi.size()));
  EXPECT_TRUE(ArraysMatch(apdq, apdq_gold, efi.size()));
}
