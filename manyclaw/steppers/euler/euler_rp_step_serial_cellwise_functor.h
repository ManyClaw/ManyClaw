#ifndef EULER_RP_STEP_SERIAL_CELLWISE_FUNCTOR_H
#define EULER_RP_STEP_SERIAL_CELLWISE_FUNCTOR_H

#include "../../riemann/euler_rp.h"
#include "../cellwise.h"

void euler_rp_step_serial_cellwise_functor( const real* q,
                                        const real* aux,
                                        const int nx,
                                        const int ny,
                                        real* amdq,
                                        real* apdq,
                                        real* wave,
                                        real* wave_speeds);

#endif // EULER_RP_STEP_SERIAL_CELLWISE_H

