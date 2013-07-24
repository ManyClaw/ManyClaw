#include <limits.h>
#include <manyclaw/manyclaw.h>
#include "gtest/gtest.h"

// Test the indexer methods
TEST(AdvectionRP, right) {
  real q_left[1] = {0.0};
  real q_right[1] = {1.0};
  real amdq[1], apdq[1], wave[1], s[1];
  advection_rp_aux_global_t aux_global;
  aux_global.u[0] = 1.0;
  aux_global.u[1] = 1.0;

  advection_rp(q_left, q_right, NULL, NULL, &aux_global, 0,
               amdq, apdq, wave, s);

  EXPECT_EQ(wave[0], 1);
  EXPECT_EQ(s[0], 1);
  EXPECT_EQ(amdq[0], 0);
  EXPECT_EQ(apdq[0], 1);
}

TEST(AdvectionRP, left) {
  real q_left[1] = {1.0};
  real q_right[1] = {0.0};
  real amdq[1], apdq[1], wave[1], s[1];
  advection_rp_aux_global_t aux_global;

  aux_global.u[0] = 1.0;
  aux_global.u[1] = 1.0;

  advection_rp(q_left, q_right, NULL, NULL, &aux_global, 0,
               amdq, apdq, wave, s);

  EXPECT_EQ(wave[0], -1);
  EXPECT_EQ(s[0], 1);
  EXPECT_EQ(amdq[0], 0);
  EXPECT_EQ(apdq[0], -1);
}
