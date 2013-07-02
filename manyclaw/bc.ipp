#ifndef __BC_IPP
#define __BC_IPP

// Zero-order Extrapolation boundary conditions
inline void set_zero_order_extrap_BCs(real* q, real* aux, const int nx, const int ny,
                               const int num_ghost, const int num_eqn)
{
    FieldIndexer fi(nx, ny, num_ghost, num_eqn);

    // Bottom edge
    for (int row = 0; row < num_ghost; ++row)
    {
        for (int col = 0; col < fi.num_col(); ++col)
        {
            for (int eqn = 0; eqn < num_eqn; ++eqn)
            {
                // q[fi.idx(row, col) + eqn] = q[fi.idx(num_ghost,col) + eqn];
                q[fi.idx(row, col) + eqn] = 1.0;
            }
        }
    }

    // Top edge
    for (int row = fi.num_row() - num_ghost; row < fi.num_row(); ++row)
    {
        for (int col = 0; col < fi.num_col(); ++col)
        {
            for (int eqn = 0; eqn < num_eqn; ++eqn)
            {
                // q[fi.idx(row, col) + eqn] = q[fi.idx(num_ghost + ny,col) + eqn];    
                q[fi.idx(row, col) + eqn] = 2.0;
            }
        }
    }

    // Left edge
    for (int row = 0; row < fi.num_row(); ++row)
    {
        for (int col = 0; col < num_ghost; ++col)
        {
            for (int eqn = 0; eqn < num_eqn; ++eqn)
            {
                // q[fi.idx(row, col) + eqn] = q[fi.idx(row, num_ghost + 1) + eqn];
                q[fi.idx(row, col) + eqn] = 3.0;
            }
        }
    }

    // Right edge
    for (int row = 0; row < fi.num_row(); ++row)
    {
        for (int col = fi.num_col() - num_ghost; col < fi.num_col(); ++col)
        {
            for (int eqn = 0; eqn < num_eqn; ++eqn)
               {
                   // q[fi.idx(row, col) + eqn] = 
                   //                 q[fi.idx(row, fi.num_col() - num_ghost - 1)];

                   q[fi.idx(row, col) + eqn] = 4.0;
               }   
        }
    }
}

// Set periodic BCs in all directions
void set_all_periodic_BCs(real* q, real* aux, const int nx, const int ny,
                               const int num_ghost, const int num_eqn)
{
    FieldIndexer fi(nx, ny, num_ghost, num_eqn);

    // Bottom edge
    for (int row = 0; row < num_ghost; ++row)
    {
        for (int col = 0; col < fi.num_col(); ++col)
        {
            for (int eqn = 0; eqn < num_eqn; ++eqn)
            {
                q[fi.idx(row, col) + eqn] = 
                           q[fi.idx(row + fi.num_row() - num_ghost, col) + eqn];
            }
        }
    }

    // Top edge
    for (int row = fi.num_row() - num_ghost; row < fi.num_row(); ++row)
    {
        for (int col = 0; col < fi.num_col(); ++col)
        {
            for (int eqn = 0; eqn < num_eqn; ++eqn)
            {
                q[fi.idx(row, col) + eqn] = 
                           q[fi.idx(row - fi.num_row() + num_ghost, col) + eqn];
            }
        }
    }

    // Left edge
    for (int row = 0; row < fi.num_row(); ++row)
    {
        for (int col = 0; col < num_ghost; ++col)
        {
            for (int eqn = 0; eqn < num_eqn; ++eqn)
            {
                q[fi.idx(row, col) + eqn] = 
                       q[fi.idx(row, col + fi.num_col() - num_ghost - 1) + eqn];
            }
        }
    }

    // Right edge
    for (int row = 0; row < fi.num_row(); ++row)
    {
        for (int col = fi.num_col() - num_ghost; col < fi.num_col(); ++col)
        {
            for (int eqn = 0; eqn < num_eqn; ++eqn)
               {
                   q[fi.idx(row, col) + eqn] = 
                                   q[fi.idx(row, col - fi.num_col() + num_ghost + 1)];
               }   
        }
    }
}

#endif