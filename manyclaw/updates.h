#ifndef __UPDATES_H
#define __UPDATES_H

#include <manyclaw/manyclaw.h>

// Upwind updater
inline
void updater_first_order_dimensional_splitting(real* q,
                                               const real* aux,
                                               const int nx,
                                               const int ny,
                                               const real* amdq,
                                               const real* apdq,
                                               const real* wave,
                                               const real* wave_speeds,
                                               const int num_ghost, 
                                               const int num_eqns,
                                               const real dtdx)
{
  int col, row, eqn;

  FieldIndexer fi(nx, ny, num_ghost, num_eqns);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqns);

  for(row = num_ghost; row <= ny + num_ghost; ++row) 
  {
    for(col = num_ghost; col <= nx + num_ghost; ++col) 
    {
      for(int eqn=0; eqn < num_eqns; ++eqn)
      {
        q[fi.left(row, col) + eqn] -= dtdx * amdq[efi.left_edge(row, col)];
        q[fi.idx(row, col) + eqn]  -= dtdx * apdq[efi.left_edge(row, col)];
        q[fi.idx(row, col) + eqn]  -= dtdx * amdq[efi.up_edge(row, col)];
        q[fi.up(row, col) + eqn]   -= dtdx * apdq[efi.up_edge(row, col)];
      }
    }
  }
}

#endif
