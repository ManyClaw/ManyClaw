#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>

#include "advection_rp_step_serial.h"
#include "advection_rp_step_tbb.h"
#include "advection_rp_step_ispc.h"
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

double benchmark(int nx, int ny, int numGhost, int numStates, int numWaves, advection_rp_step_t rp_stepper)
{
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

  if (argc == 3)
  {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
  }

  printf("Riemann Solve on %dx%d grid\n", nx, ny);

  const char * advection_rp_stepper_names[] =
    {
      "serial",
      "TBB"
    };

  advection_rp_step_t advection_rp_steppers[] =
    {
      advection_rp_step_serial,
      advection_rp_step_tbb
      // TODO add other advection_rp_step functions here
    };

  size_t num_rp_kernels = sizeof(advection_rp_steppers)/sizeof(advection_rp_step_t);

  for (size_t i = 0; i < num_rp_kernels; i++)
  {
    std::cout << "Benchmarking " << advection_rp_stepper_names[i] << " Riemann kernel...\n";
    std::cout << "  finished in " << 1e3 * benchmark(nx, ny, 2, 3, 2, advection_rp_steppers[i]) << " ms\n";
  }
  
  return 0;
}

