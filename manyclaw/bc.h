#ifndef BC_H
#define BC_H

typedef double real;

// Zero-order Extrapolation boundary conditions
inline
void set_zero_order_extrap_BCs(real* q, real* aux, const int nx, const int ny,
                               const int num_ghost, const int num_eqn)
{
    int idx_gc, idx_val;
    int num_in_row = (2 * num_ghost + nx) * num_eqn;
    // int num_in_col = (2 * num_ghost + ny) * num_eqn;

    // Bottom edge
    // for (int i=0; i < nx; ++i)
    //   for (int j=0; j < num_ghost; ++j)
    //     for (int m=0; m < num_eqn; ++m)
    //       q[m + i * num_eqn + j * num_eqn * (2 * num_ghost + nx)] = 100.0;
    //                       // q[m + i * num_eqn + num_ghost * (2 * num_ghost + nx)];
    for (int col = num_ghost; col < num_ghost + nx; ++col)
    {
      for (int row = 0; row < num_ghost; ++row)
      {
        idx_gc = num_in_row * row + col * num_eqn;
        idx_val = idx_gc + num_in_row * (num_ghost - row);
        for (int eqn = 0; eqn < num_eqn; ++eqn)
        {
          q[idx_gc + eqn] = q[idx_val + eqn];
          q[idx_gc + eqn] = 200.0;
        }
      }
    }

    // Top edge
    // for (int i = 0; i < nx; ++i)
    //   for (int j = 0; j < num_ghost; ++j)
    //     for (int m = 0; m < num_eqn; ++m)
    //       q[m + i * num_eqn + (j + num_ghost + ny + 1) * (2 * num_ghost + nx)] = 200.0;
              // q[m + i * num_eqn + (num_ghost + ny) * (2 * num_ghost + nx)];
    // for (int col = num_ghost; col < num_ghost + nx; ++col)
    // {
    //   for (int row = 0; row < 2 * num_ghost + ny; ++row)
    //   {
    //     idx_gc = num_in_row * row + col * num_eqn;
    //     idx_val = idx_gc - num_in_row * (num_ghost - row);
    //     for (int eqn = 0; eqn < num_eqn; ++eqn)
    //     {
    //       q[idx_gc + eqn] = q[idx_val + eqn];
    //     }
    //   }
    // }

    // Left and right edges
    // for (int i = 0; i < num_ghost; ++i)
    // {
    //   for (int j = 0; j < ny; ++j)
    //   {
    //     for (int m = 0; m < num_eqn; ++m)
    //     {
    //       q[m + i * num_eqn + j * (2 * num_ghost + nx)] = 300.0;
    //                 // q[m + num_ghost * num_eqn + j * (2 * num_ghost + nx)];
    //       q[m + (i + num_ghost + nx) + j * (2 * num_ghost + nx)] = 400.0;
    //                 // q[m + (num_ghost + nx - 1) + j * (2 * num_ghost + nx)];
    //     }
    //   }
    // }

    // Left
    for (int row = num_ghost; row < num_ghost + ny; ++row)
    {
      for (int col = 0; col < num_ghost; ++col)
      {
        idx_gc = num_in_row * row + col * num_eqn;
        idx_val = idx_gc + (num_ghost - col);
        for (int eqn = 0; eqn < num_eqn; ++eqn)
        {
          q[idx_gc + eqn] = q[idx_val + eqn];
          q[idx_gc + eqn] = 100.0;
        }
      }
    }
      
}

// Set periodic BCs in all directions
inline
void set_all_periodic_BCs(real* q, real* aux, const int nx, const int ny,
                               const int num_ghost, const int num_states)
{

}

#endif
