#include "euler_rp_step_serial_cellwise_functor.h"

void euler_rp_step_serial_cellwise_functor(const real* q,  const real* aux,
                                        const int nx, const  int ny,
                                        real* amdq, real* apdq, real* wave,
                                        real* wave_speeds)
{
  RiemannCellwiseStepper()(euler_rp, euler_rp_grid_params,
                           q, aux, &euler_rp_aux_global, nx, ny,
                           amdq, apdq, wave, wave_speeds);
}

