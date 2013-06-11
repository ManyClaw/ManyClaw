#ifndef __BC_H
#define __BC_H

typedef double real;

// Zero-order Extrapolation boundary conditions
void set_zero_order_extrap_BCs(real* q, real* aux, const int nx, const int ny,
                               const int num_ghost, const int num_eqn)
{
    // Bottom edge
    for (int i=0; i < nx; ++i)
      for (int j=0; j < num_ghost; ++j)
        for (int m=0; m < num_eqn; ++m)
          q[m + i * num_eqn + j * (2 * num_ghost + nx)] = 100.0;
                          // q[m + i * num_eqn + num_ghost * (2 * num_ghost + nx)];
    // Top edge
    for (int i = 0; i < nx; ++i)
      for (int j = 0; j < num_ghost; ++j)
        for (int m = 0; m < num_eqn; ++m)
          q[m + i * num_eqn + (j + num_ghost + ny + 1) * (2 * num_ghost + nx)] = 200.0;
              // q[m + i * num_eqn + (num_ghost + ny) * (2 * num_ghost + nx)];

    // Left and right edges
    for (int i = 0; i < num_ghost; ++i)
    {
      for (int j = 0; j < ny; ++j)
      {
        for (int m = 0; m < num_eqn; ++m)
        {
          q[m + i * num_eqn + j * (2 * num_ghost + nx)] = 300.0;
                    // q[m + num_ghost * num_eqn + j * (2 * num_ghost + nx)];
          q[m + (i + num_ghost + nx) + j * (2 * num_ghost + nx)] = 400.0;
                    // q[m + (num_ghost + nx - 1) + j * (2 * num_ghost + nx)];
        }
      }
    }
      
}

// Set periodic BCs in all directions
void set_all_periodic_BCs(real* q, real* aux, const int nx, const int ny,
                               const int num_ghost, const int num_states)
{

}

#endif