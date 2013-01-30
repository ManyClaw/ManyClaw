#ifndef COMMON_H
#define COMMON_H

#include "../many_claw.h"


struct rp_grid_params
{
  int num_ghost;
  int num_states;
  int num_aux;
  int num_waves;
};

typedef void (*rp_t)(const real* q_left, const real* q_right,
                  const real* aux_left, const real* aux_right,
                  const void* aux_global,
                     real* amdq, real* apdq, real* wave, real* s);

typedef void (*rp_step_t)(const real* q, const real* aux,
                          const int nx, const int ny, real* amdq,  real* apdq,
                          real* wave, real* wave_speeds);

typedef void (*updater_t)(real* q, const real* aux, const int nx, const int ny,
                                   const real* amdq, const real* apdq,
                                   const real* wave, const real* wave_speeds,
                                   const rp_grid_params grid_params);

struct abs_diff
{
  real operator()(real a, real b) const
  {
    return std::abs(a - b);
  }
};

template <typename Vector>
void randomize_vector(Vector& v)
{
  // generate numbers in [0,1) deterministically
  real seed = 0.123456789;

  for(size_t i = 0; i < v.size(); i++)
  {
    seed = (314159.26535 * seed);
    seed = seed - floor(seed);
    v[i] = seed;
  }
}

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

void compare_steppers(int nx, int ny, rp_grid_params params,
                      rp_step_t rp_stepper1, rp_step_t rp_stepper2);
double benchmark_stepper(int nx, int ny, rp_grid_params params,
                         rp_step_t rp_stepper);
double compare_updates(int nx, int ny, rp_grid_params params, 
                         rp_step_t rp_stepper, updater_t updater);

#endif // COMMON_H
