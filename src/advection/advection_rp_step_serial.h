#ifndef ADVECTION_RP_STEP_SERIAL_H
#define ADVECTION_RP_STEP_SERIAL_H

#include "advection_rp.h"

void advection_rp_step_serial(real* q,     real* aux,
                              int numGhost, int numStates, int numWaves, int nx, int ny, // inputs
                              real* amdq,  real* apdq,
                              real* wave,  real* wave_speeds);

#endif // ADVECTION_RP_STEP_SERIAL_H

