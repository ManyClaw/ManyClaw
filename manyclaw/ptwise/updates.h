#ifndef __UPDATES_H
#define __UPDATES_H

typedef double real;

// Upwind updater
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
  int col, row, idx_left, idx_center, idx_up, idx_out_x, idx_out_y;

#pragma omp parallel for schedule(runtime)
  for(row = num_ghost; row <= ny + num_ghost; ++row) 
  {
    for(col = num_ghost; col <= nx + num_ghost; ++col) 
    {
      idx_left = col + row*(nx + 2*num_ghost) - 1;
      idx_up = col + (row - 1)*(nx + 2*num_ghost);
      idx_center = idx_left + 1;
      idx_out_x = (col - num_ghost) + (row - num_ghost) * (nx + 1);
      idx_out_y = idx_out_x + ((nx + 1)*(ny + 1));

      printf("center, up, left: %d, %d, %d\n", idx_center, idx_up, idx_left);
      printf("  out_x, out_y: %d, %d\n", idx_out_x, idx_out_y);

      for(int eqn=0; eqn < num_eqns; ++eqn)
      {
        q[idx_left*num_eqns + eqn]   -= dtdx * amdq[idx_out_x*num_eqns + eqn];
        q[idx_up*num_eqns + eqn]     -= dtdx * apdq[idx_out_y*num_eqns + eqn];
        q[idx_center*num_eqns + eqn] -= dtdx * apdq[idx_out_x*num_eqns + eqn];
        q[idx_center*num_eqns + eqn] -= dtdx * amdq[idx_out_y*num_eqns + eqn];
      }
    }
  }
}

#endif
