#include <stdlib.h>
#include <stdio.h>

// TODO add other step headers here
#include "advection_rp_step_serial.h"

void benchmark(int nx, int ny, int numGhost, int numStates, int numWaves)
{
  real* q     = (real*) malloc(sizeof(real)*(nx+numGhost*2)*(ny+numGhost*2)*numStates);
  real* aux   = (real*) malloc(sizeof(real)*(nx+numGhost*2)*(ny+numGhost*2)*2);
  real* amdq  = (real*) malloc(sizeof(real)*(nx+numGhost*2)*(ny+numGhost*2)*numStates);
  real* apdq  = (real*) malloc(sizeof(real)*(nx+numGhost*2)*(ny+numGhost*2)*numStates);
  real* waves = (real*) malloc(sizeof(real)*(nx+numGhost*2)*(ny+numGhost*2)*numStates*numWaves);
  real* wave_speeds = (real*) malloc(sizeof(real)*(nx+numGhost*2)*(ny+numGhost*2)*4*numStates);

  // TODO start timer

  advection_rp_step_serial(q, aux, numGhost, numStates, numWaves, nx, ny,
                           amdq, apdq, waves, wave_speeds);
  
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

  benchmark(nx, ny, 2, 3, 2);
  
  return 0;
}

