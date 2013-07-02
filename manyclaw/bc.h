#ifndef __BC_H
#define __BC_H

#include "../common/common.h"

typedef void (*set_bc_t)(real* q, real* aux, const int nx, const int ny,
                               const int num_ghost, const int num_eqn);

#include "bc.ipp"

#endif