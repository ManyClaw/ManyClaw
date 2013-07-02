#ifndef __LIMITERS_H
#define __LIMITERS_H

#include <algorithm>
#include <cmath>

typedef double real;

inline void limit_waves(int num_eqn, int num_waves, int num_ghost, 
                        int num_cells, real* wave, real* s, 
                        real (*limiter_func)(real))
{
    real wave_norm_squared = 0;
    real left_dot_product = 0;
    real right_dot_product = 0;
    real r, phi;
    real const epsilon = 1e-90;

    for (int wave_index = 0; wave_index < num_waves; wave_index++)
    {
        right_dot_product = 0.0;
        for (int i = num_ghost - 1; i < num_cells + num_ghost + 1; i++)
        {

            wave_norm_squared = 0.0;
            left_dot_product = right_dot_product;
            right_dot_product = 0.0;
            for (int eqn_index = 0; eqn_index < num_eqn; eqn_index++)
            {
                wave_norm_squared += pow(wave[eqn_index 
                                            + wave_index * num_eqn 
                                            + i * num_waves * num_eqn],2);
                right_dot_product += wave[eqn_index 
                                        + wave_index * num_eqn 
                                        + i * num_waves * num_eqn] 
                                   * wave[eqn_index 
                                        + wave_index * num_eqn 
                                        + (i+1) * num_waves * num_eqn];
            }
            
            // Here to just calculate the left_dot_product, skip the iteration
            if (i == 0) continue;

            // If wave norm is too small, skip this cell
            if (wave_norm_squared < epsilon) continue;

            // Evaluate limiter
            if (s[wave_index * num_eqn + i * num_waves * num_eqn] > 0.0)
                r = left_dot_product / wave_norm_squared;
            else
                r = right_dot_product / wave_norm_squared;
            phi = limiter_func(r);

            // Apply limiter to this wave
            for (int eqn_index = 0; eqn_index < num_eqn; eqn_index++)
                wave[eqn_index 
                   + wave_index * num_eqn 
                   + i * num_waves * num_eqn] *= phi;
                
        }
    }
}

/* In classic Clawpack the limiters are labeled by number via
    no limiting = 0
    std::minmod = 1
    superbee = 2
    vanleer = 3
    mc = 4
    beam-warstd::ming = 5
*/

// MinMod
inline real minmod_limiter(real r)
{
    return std::max(0.0,std::min(1.0,r));
}

// Superbee
inline real superbee_limiter(real r)
{
    return std::max(0.0,std::max(std::min(1.0,2.0*r),std::min(2.0,r)));
}

// Van Leer
inline real vanleer_limiter(real r)
{
    return (r + std::abs(r)) / (1.0 + std::abs(r));
}

// Monotonized Center
inline real mc_limiter(real r)
{
    return std::max(0.0,std::min((1.0 + r) / 2.0,std::min(2.0,2.0*r)));
}

// Beam-Warstd::ming
inline real beamwarming_limiter(real r)
{
    return r;
}

#endif
