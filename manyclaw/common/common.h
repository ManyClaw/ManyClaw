#ifndef __COMMON_H
#define __COMMON_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "data_structures.h"

struct abs_diff
{
  real operator()(real a, real b) const
  {
    return std::abs(a - b);
  }
};

// maximum absolute difference between components
template <typename Vector>
real max_error(const Vector& v1, const Vector& v2)
{
  return std::inner_product
    (v1.begin(), v1.end(),
     v2.begin(),
     real(0),
     std::max<real>,
     abs_diff());
}

int test_function(int nx);

void compare_grid_evals(int nx, int ny, rp_grid_params_t params,
                        rp_grid_eval_t rp_grid_eval1, rp_grid_eval_t rp_grid_eval2,
                        void* aux_global);
double benchmark_grid_eval(int nx, int ny, rp_grid_params_t params,
                           rp_grid_eval_t rp_grid_eval,
                           void* aux_global);
double compare_updates(int nx, int ny, rp_grid_params_t params, 
                       rp_grid_eval_t rp_grid_eval, updater_t updater);

#endif // COMMON_H
