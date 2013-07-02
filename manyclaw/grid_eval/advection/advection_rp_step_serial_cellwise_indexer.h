#ifndef ADVECTION_RP_STEP_SERIAL_CELLWISE_INDEXER_H
#define ADVECTION_RP_STEP_SERIAL_CELLWISE_INDEXER_H

#include "../../common/data_structures.h"
#include "../../riemann/advection_rp.h"

void advection_rp_step_serial_cellwise_indexer( const real* q,
                                        const real* aux,
                                        const int nx,
                                        const int ny,
                                        real* amdq,
                                        real* apdq,
                                        real* wave,
                                        real* wave_speeds);

#endif // ADVECTION_RP_STEP_SERIAL_CELLWISE_INDEXER_H

