#ifndef __BC_IPP
#define __BC_IPP

// Zero-order Extrapolation boundary conditions
inline void set_zero_order_extrap_BCs(real* q, real* aux, const unsigned nx, const unsigned ny,
                               const unsigned num_ghost, const unsigned num_eqn)
{
    FieldIndexer fi(nx, ny, num_ghost, num_eqn);

    // Bottom edge
    for (unsigned row = 0; row < num_ghost; ++row)
    {
        for (unsigned col = 0; col < fi.num_col(); ++col)
        {
            for (unsigned eqn = 0; eqn < num_eqn; ++eqn)
            {
                q[fi.idx(row, col) + eqn] = q[fi.idx(num_ghost,col) + eqn];
            }
        }
    }

    // Top edge
    for (unsigned row = fi.num_row() - num_ghost; row < fi.num_row(); ++row)
    {
        for (unsigned col = 0; col < fi.num_col(); ++col)
        {
            for (unsigned eqn = 0; eqn < num_eqn; ++eqn)
            {
                q[fi.idx(row, col) + eqn] = q[fi.idx(num_ghost + ny,col) + eqn];    
            }
        }
    }

    // Left edge
    for (unsigned row = 0; row < fi.num_row(); ++row)
    {
        for (unsigned col = 0; col < num_ghost; ++col)
        {
            for (unsigned eqn = 0; eqn < num_eqn; ++eqn)
            {
                q[fi.idx(row, col) + eqn] = q[fi.idx(row, num_ghost + 1) + eqn];
            }
        }
    }

    // Right edge
    for (unsigned row = 0; row < fi.num_row(); ++row)
    {
        for (unsigned col = fi.num_col() - num_ghost; col < fi.num_col(); ++col)
        {
            for (unsigned eqn = 0; eqn < num_eqn; ++eqn)
               {
                   q[fi.idx(row, col) + eqn] = 
                                   q[fi.idx(row, fi.num_col() - num_ghost - 1)];
               }   
        }
    }
}

// Set periodic BCs in all directions
inline void set_all_periodic_BCs(real* q, real* aux, const unsigned nx, const unsigned ny,
                               const unsigned num_ghost, const unsigned num_eqn)
{
    FieldIndexer fi(nx, ny, num_ghost, num_eqn);

    // Bottom edge
    for (unsigned row = 0; row < num_ghost; ++row)
    {
        for (unsigned col = 0; col < fi.num_col(); ++col)
        {
            for (unsigned eqn = 0; eqn < num_eqn; ++eqn)
            {
                q[fi.idx(row, col) + eqn] = 
                           q[fi.idx(row + fi.num_row() - num_ghost, col) + eqn];
            }
        }
    }

    // Top edge
    for (unsigned row = fi.num_row() - num_ghost; row < fi.num_row(); ++row)
    {
        for (unsigned col = 0; col < fi.num_col(); ++col)
        {
            for (unsigned eqn = 0; eqn < num_eqn; ++eqn)
            {
                q[fi.idx(row, col) + eqn] = 
                           q[fi.idx(row - fi.num_row() + num_ghost, col) + eqn];
            }
        }
    }

    // Left edge
    for (unsigned row = 0; row < fi.num_row(); ++row)
    {
        for (unsigned col = 0; col < num_ghost; ++col)
        {
            for (unsigned eqn = 0; eqn < num_eqn; ++eqn)
            {
                q[fi.idx(row, col) + eqn] = 
                       q[fi.idx(row, col + fi.num_col() - num_ghost - 1) + eqn];
            }
        }
    }

    // Right edge
    for (unsigned row = 0; row < fi.num_row(); ++row)
    {
        for (unsigned col = fi.num_col() - num_ghost; col < fi.num_col(); ++col)
        {
            for (unsigned eqn = 0; eqn < num_eqn; ++eqn)
               {
                   q[fi.idx(row, col) + eqn] = 
                                   q[fi.idx(row, col - fi.num_col() + num_ghost + 1)];
               }   
        }
    }
}

#endif
