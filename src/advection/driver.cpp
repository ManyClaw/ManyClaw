#include <stdlib.h>
#include <stdio.h>
#include <iostream>

// TODO add other step headers here
#include "advection_rp_step_serial.h"
#include "advection_rp_step_tbb.h"
#include "advection_rp_step_ispc.h"

void benchmark(int nx, int ny, int numGhost, int numStates, int numWaves, advection_rp_step_t rp_stepper)
{
  real* q     = (real*) malloc(sizeof(real)*(nx+numGhost*2)*(ny+numGhost*2)*numStates);
  real* aux   = (real*) malloc(sizeof(real)*(nx+numGhost*2)*(ny+numGhost*2)*2);
  real* amdq  = (real*) malloc(sizeof(real)*(nx+numGhost*2)*(ny+numGhost*2)*numStates);
  real* apdq  = (real*) malloc(sizeof(real)*(nx+numGhost*2)*(ny+numGhost*2)*numStates);
  real* waves = (real*) malloc(sizeof(real)*(nx+numGhost*2)*(ny+numGhost*2)*numStates*numWaves);
  real* wave_speeds = (real*) malloc(sizeof(real)*(nx+numGhost*2)*(ny+numGhost*2)*4*numStates);

  // TODO start timer

  rp_stepper(q, aux, numGhost, numStates, numWaves, nx, ny, amdq, apdq, waves, wave_speeds);
  
  // TODO stop timer
  
  free(q);
  free(aux);
  free(waves);
  free(amdq);
  free(apdq);
  free(wave_speeds);
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
    std::cout << "Benchmarking " << advection_rp_stepper_names[i] << " Riemann kernel\n";
    benchmark(nx, ny, 2, 3, 2, advection_rp_steppers[i]);
  }
  
  return 0;
}

