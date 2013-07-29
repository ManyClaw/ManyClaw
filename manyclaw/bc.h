#ifndef __BC_H
#define __BC_H

#include "common/common.h"

typedef void (*set_bc_t)(real* q, real* aux, const unsigned nx, const unsigned ny,
                               const unsigned num_ghost, const unsigned num_eqn);

#include "bc.ipp"

#endif
