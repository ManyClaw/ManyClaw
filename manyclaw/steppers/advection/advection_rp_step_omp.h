#ifndef ADVECTION_RP_STEP_OMP_H
#define ADVECTION_RP_STEP_OMP_H

#include "../../riemann/advection_rp.h"

void advection_rp_step_omp(const real* q,
                           const real* aux,
                           const int nx,
                           const int ny,
                           real* amdq,
                           real* apdq,
                           real* wave,
                           real* wave_speeds);

#endif // ADVECTION_RP_STEP_OMP_H
