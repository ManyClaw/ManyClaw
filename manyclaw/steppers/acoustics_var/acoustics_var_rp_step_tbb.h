#ifndef ACOUSTICS_VAR_RP_STEP_TBB_H
#define ACOUSTICS_VAR_RP_STEP_TBB_H

#include "../../riemann/acoustics_var_rp.h"

void acoustics_var_rp_step_tbb(const real* q,
                           const real* aux,
                           const int nx,
                           const int ny,
                           real* amdq,
                           real* apdq,
                           real* wave,
                           real* wave_speeds);

#endif // ACOUSTICS_VAR_RP_STEP_TBB_H

