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
                                               const int num_eqn,
                                               const real dtdx)
{
  FieldIndexer fi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn);

  for(int row = num_ghost; row <= ny + num_ghost; ++row) 
  {
    for(int col = num_ghost; col <= nx + num_ghost; ++col) 
    {
      for(int eqn=0; eqn < num_eqn; ++eqn)
      {
        q[fi.left(row, col) + eqn] -= dtdx * amdq[efi.left_edge(row, col) + eqn];
        q[fi.idx(row, col) + eqn]  -= dtdx * apdq[efi.left_edge(row, col) + eqn];
        q[fi.down(row, col) + eqn] -= dtdx * amdq[efi.down_edge(row, col) + eqn];
        q[fi.idx(row, col) + eqn]  -= dtdx * apdq[efi.down_edge(row, col) + eqn];
      }
    }
  }
}

#endif
