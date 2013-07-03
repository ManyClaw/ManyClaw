#ifndef __MANYCLAW_TEST_UTILS_H
#define __MANYCLAW_TEST_UTILS_H

#include "gtest/gtest.h"

#include <manyclaw/manyclaw.h>

#include <limits.h>


::testing::AssertionResult ArraysMatch(const real *expected,
                                       const real *actual, int size);

::testing::AssertionResult GridEvalsMatch(int nx, int ny, rp_grid_params params,
					  rp_grid_eval_t rp_grid_eval1, 
					  rp_grid_eval_t rp_grid_eval2);

double compare_updates(int nx, int ny, rp_grid_params params, 
                       rp_grid_eval_t rp_grid_eval, 
                       updater_t updater);

#endif // __MANYCLAW_TEST_UTILS_H
