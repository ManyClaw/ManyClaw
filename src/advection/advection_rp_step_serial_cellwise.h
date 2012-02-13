#ifndef ADVECTION_RP_STEP_SERIAL_CELLWISE_H
#define ADVECTION_RP_STEP_SERIAL_CELLWISE_H

#include "advection_rp.h"

void advection_rp_step_serial_cellwise( real* q,  real* aux,  int num_ghost,
   int num_states,  int num_waves,  int nx,  int ny,
  real* amdq, real* apdq, real* wave, real* wave_speeds);

#endif // ADVECTION_RP_STEP_SERIAL_CELLWISE_H

