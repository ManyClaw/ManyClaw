#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

#include "../ptwise/riemann/advection_rp.h"
#include "advection_rp_step_serial.h"
//#include "advection_rp_step_serial_tiled.h"
//#include "advection_rp_step_serial_cellwise.h"
//#include "advection_rp_step_tbb.h"
//#include "advection_rp_step_ispc.h"
//#include "advection_rp_step_omp.h"
// TODO add other step headers here

#include "timer.h"

template <typename Vector>
void randomize(Vector& v)
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


void compare_kernels(int nx, int ny, advection_rp_step_t rp_stepper1,
                     advection_rp_step_t rp_stepper2)
{
  const int num_aux = 2;
  const int dim = 2;
  const int num_ghost = advection_rp_params.num_ghost;
  const int num_states = advection_rp_params.num_states;
  const int num_waves  = advection_rp_params.num_waves;

  // TODO make this a function of the encapsulated states
  // Cell centered values
  std::vector<real> q           ((nx+num_ghost*2)*(ny+num_ghost*2)*num_states);
  std::vector<real> aux         ((nx+num_ghost*2)*(ny+num_ghost*2)*num_aux);

  randomize(q);
  randomize(aux);

  // Outputs on interfaces
  std::vector<real> amdq1        ((nx+1)*(ny+1)*num_states*dim);
  std::vector<real> apdq1        ((nx+1)*(ny+1)*num_states*dim);
  std::vector<real> waves1       ((nx+1)*(ny+1)*num_states*num_waves*dim);
  std::vector<real> wave_speeds1 ((nx+1)*(ny+1)*4*num_states*dim);

  std::vector<real> amdq2        ((nx+1)*(ny+1)*num_states*dim);
  std::vector<real> apdq2        ((nx+1)*(ny+1)*num_states*dim);
  std::vector<real> waves2       ((nx+1)*(ny+1)*num_states*num_waves*dim);
  std::vector<real> wave_speeds2 ((nx+1)*(ny+1)*4*num_states*dim);

  rp_stepper1(&q[0], &aux[0], nx, ny,
              &amdq1[0], &apdq1[0], &waves1[0], &wave_speeds1[0]);
  rp_stepper2(&q[0], &aux[0], nx, ny,
              &amdq2[0], &apdq2[0], &waves2[0], &wave_speeds2[0]);

  std::cout << "  Comparing solution to reference solution...\n";
  std::cout << "    amdq        " << max_error(amdq1, amdq2) << "\n";
  std::cout << "    apdq        " << max_error(apdq1, apdq2) << "\n";
  std::cout << "    waves       " << max_error(waves1, waves2) << "\n";
  std::cout << "    wave_speeds " << max_error(wave_speeds1, wave_speeds2) << "\n";
}

double benchmark(int nx, int ny, advection_rp_step_t rp_stepper)
{
  const int num_aux = 2;
  const int dim = 2;
  const int num_ghost = advection_rp_params.num_ghost;
  const int num_states = advection_rp_params.num_states;
  const int num_waves  = advection_rp_params.num_waves;

  // TODO encapsulate all state into a single structure
  std::vector<real> q           ((nx+num_ghost*2)*(ny+num_ghost*2)*num_states);
  std::vector<real> aux         ((nx+num_ghost*2)*(ny+num_ghost*2)*num_aux);
  std::vector<real> amdq        (((nx+1)*(ny+1)*num_states)*dim);
  std::vector<real> apdq        (((nx+1)*(ny+1)*num_states)*dim);
  std::vector<real> waves       ((nx+1)*(ny+1)*num_states*num_waves*dim);
  std::vector<real> wave_speeds ((nx+1)*(ny+1)*4*num_states*dim);

  randomize(q);
  randomize(aux);

  timer t;

  rp_stepper(&q[0], &aux[0], nx, ny, &amdq[0], &apdq[0], &waves[0], &wave_speeds[0]);

  double time = t.elapsed();

  return time;
}

int main(int argc, char ** argv)
{
  int nx = 100, ny = 100;

  if (argc == 3)
  {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
  }

  printf("Riemann Solve on %dx%d grid\n", nx, ny);

  const char * advection_rp_stepper_names[] =
    {
      "serial",
      "serial_tiled",
      "serial_cellwise",
      "TBB",
      "omp"
    };

  advection_rp_step_t advection_rp_steppers[] =
    {
      advection_rp_step_serial,
      //advection_rp_step_serial_tiled,
      //advection_rp_step_serial_cellwise,
      //advection_rp_step_tbb,
      //advection_rp_step_omp
      // TODO add other advection_rp_step functions here
    };

  size_t num_rp_kernels = sizeof(advection_rp_steppers)/sizeof(advection_rp_step_t);

  for (size_t i = 0; i < num_rp_kernels; i++)
  {
    //std::cout << "Testing " << advection_rp_stepper_names[i] << " Riemann kernel...\n";
    //compare_kernels(nx, ny, advection_rp_step_serial, advection_rp_steppers[i]);
    std::cout << "  Benchmark finished in " << 1e3 * benchmark(nx, ny, advection_rp_steppers[i]) << " ms\n";
  }

  return 0;
}

