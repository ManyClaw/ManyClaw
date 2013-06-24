#include <limits.h>
#include <manyclaw/manyclaw.h>
#include "gtest/gtest.h"

// Test the indexer methods
TEST(AdvectionRP, base) {
  real q_left[1] = {0.0};
  real q_right[1] = {1.0};
  real amdq[1], apdq[1], wave[1], s[1];

  advection_rp(q_left, q_right, NULL, NULL, &advection_rp_aux_global,
               amdq, apdq, wave, s);
  

  EXPECT_EQ(wave[0], 1);
  EXPECT_EQ(s[0], 1);
  EXPECT_EQ(amdq[0], 0);
  EXPECT_EQ(apdq[0], 1);
}
