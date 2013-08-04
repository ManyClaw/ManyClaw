#ifndef __MANYCLAW_TEST_UTILS_H
#define __MANYCLAW_TEST_UTILS_H

#include "gtest/gtest.h"

#include <manyclaw/manyclaw.h>

#include <limits.h>


::testing::AssertionResult ArraysMatch(const real *expected,
                                       const real *actual, int size);

::testing::AssertionResult GridEvalsMatch(int nx, int ny, rp_grid_params_t params,
                                          rp_grid_eval_t rp_grid_eval_1, const char* rp_grid_eval_name_1, 
                                          rp_grid_eval_t rp_grid_eval_2, const char* rp_grid_eval_name_2, 
                                          void* aux_global);

double compare_updates(int nx, int ny, rp_grid_params_t params, 
                       rp_grid_eval_t rp_grid_eval, 
                       updater_t updater, void* aux_global);

#endif // __MANYCLAW_TEST_UTILS_H
