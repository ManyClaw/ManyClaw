#ifndef EULER_RP_STEP_SERIAL_H
#define EULER_RP_STEP_SERIAL_H

#include "../../riemann/euler_rp.h"

void euler_rp_step_serial( const real* q,
                               const real* aux,
                               const int nx,
                               const int ny,
                               real* amdq,
                               real* apdq,
                               real* wave,
                               real* wave_speeds);
#endif // EULER_RP_STEP_SERIAL_H

