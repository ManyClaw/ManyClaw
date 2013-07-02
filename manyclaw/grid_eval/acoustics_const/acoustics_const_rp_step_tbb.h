#ifndef ACOUSTICS_CONST_RP_STEP_TBB_H
#define ACOUSTICS_CONST_RP_STEP_TBB_H

#include "../../riemann/acoustics_const_rp.h"

void acoustics_const_rp_step_tbb(const real* q,
                           const real* aux,
                           const int nx,
                           const int ny,
                           real* amdq,
                           real* apdq,
                           real* wave,
                           real* wave_speeds);

#endif // ACOUSTICS_CONST_RP_STEP_TBB_H

