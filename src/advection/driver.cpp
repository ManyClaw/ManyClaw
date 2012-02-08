#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

#include "advection_rp_step_serial.h"
#include "advection_rp_step_tbb.h"
#include "advection_rp_step_ispc.h"
#include "advection_rp_step_omp.h"
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


void compare_kernels(int nx, int ny, int numGhost, int numStates, int numWaves,
                     advection_rp_step_t rp_stepper1,
                     advection_rp_step_t rp_stepper2)
{
  // TODO make this a function of the encapsulated states
  std::vector<real> q           ((nx+numGhost*2)*(ny+numGhost*2)*numStates);
  std::vector<real> aux         ((nx+numGhost*2)*(ny+numGhost*2)*2);

  randomize(q);
  randomize(aux);

  std::vector<real> amdq1        ((nx+numGhost*2)*(ny+numGhost*2)*numStates);
  std::vector<real> apdq1        ((nx+numGhost*2)*(ny+numGhost*2)*numStates);
  std::vector<real> waves1       ((nx+numGhost*2)*(ny+numGhost*2)*numStates*numWaves);
  std::vector<real> wave_speeds1 ((nx+numGhost*2)*(ny+numGhost*2)*4*numStates);
  
  std::vector<real> amdq2        ((nx+numGhost*2)*(ny+numGhost*2)*numStates);
  std::vector<real> apdq2        ((nx+numGhost*2)*(ny+numGhost*2)*numStates);
  std::vector<real> waves2       ((nx+numGhost*2)*(ny+numGhost*2)*numStates*numWaves);
  std::vector<real> wave_speeds2 ((nx+numGhost*2)*(ny+numGhost*2)*4*numStates);

  rp_stepper1(&q[0], &aux[0], numGhost, numStates, numWaves, nx, ny,
              &amdq1[0], &apdq1[0], &waves1[0], &wave_speeds1[0]);
  rp_stepper2(&q[0], &aux[0], numGhost, numStates, numWaves, nx, ny,
              &amdq2[0], &apdq2[0], &waves2[0], &wave_speeds2[0]);

  std::cout << "  Comparing solution to reference solution...\n";
  std::cout << "    amdq        " << max_error(amdq1, amdq2) << "\n";
  std::cout << "    apdq        " << max_error(apdq1, apdq2) << "\n";
  std::cout << "    waves       " << max_error(waves1, waves2) << "\n";
  std::cout << "    wave_speeds " << max_error(wave_speeds1, wave_speeds2) << "\n";
}

double benchmark(int nx, int ny, int numGhost, int numStates, int numWaves, advection_rp_step_t rp_stepper)
{
  // TODO encapsulate all state into a single structure
  std::vector<real> q           ((nx+numGhost*2)*(ny+numGhost*2)*numStates);
  std::vector<real> aux         ((nx+numGhost*2)*(ny+numGhost*2)*2);
  std::vector<real> amdq        ((nx+numGhost*2)*(ny+numGhost*2)*numStates);
  std::vector<real> apdq        ((nx+numGhost*2)*(ny+numGhost*2)*numStates);
  std::vector<real> waves       ((nx+numGhost*2)*(ny+numGhost*2)*numStates*numWaves);
  std::vector<real> wave_speeds ((nx+numGhost*2)*(ny+numGhost*2)*4*numStates);

  randomize(q);
  randomize(aux);

  timer t;

  rp_stepper(&q[0], &aux[0], numGhost, numStates, numWaves, nx, ny, &amdq[0], &apdq[0], &waves[0], &wave_speeds[0]);

  double time = t.elapsed();

  return time;
}

int main(int argc, char ** argv)
{
  int nx = 100, ny = 100;
  int numGhost  = 2;
  int numStates = 3;
  int numWaves  = 2;

  if (argc == 3)
  {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
  }

  printf("Riemann Solve on %dx%d grid\n", nx, ny);

  const char * advection_rp_stepper_names[] =
    {
      "serial",
      //"TBB",
      "omp"
    };

  advection_rp_step_t advection_rp_steppers[] =
    {
      advection_rp_step_serial,
      //advection_rp_step_tbb,
      advection_rp_step_omp
      // TODO add other advection_rp_step functions here
    };

  size_t num_rp_kernels = sizeof(advection_rp_steppers)/sizeof(advection_rp_step_t);

  for (size_t i = 0; i < num_rp_kernels; i++)
  {
    std::cout << "Testing " << advection_rp_stepper_names[i] << " Riemann kernel...\n";
    compare_kernels(nx, ny, numGhost, numStates, numWaves, advection_rp_step_serial, advection_rp_steppers[i]);
    std::cout << "  Benchmark finished in " << 1e3 * benchmark(nx, ny, numGhost, numStates, numWaves, advection_rp_steppers[i]) << " ms\n";
  }

  return 0;
}

