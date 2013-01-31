#ifndef __LIMITERS_H
#define __LIMITERS_H

#include <algorithm>
#include <cmath>

typedef double real;

// Limit wave argument based on wave_left, wave_right, and s
inline void limiter(const int num_eqn, const int num_waves,
                    const real* wave_left, const real* wave_right,
                    real (*limiter_func)(real),
                    real* wave, real* s)
{
    int eqn_index, wave_index;
    real wave_norm_squared, left_dot_product, right_dot_product;
    real r;
    
    for (wave_index=0; wave_index < num_waves; wave_index++)
    {
        for (eqn_index=0; eqn_index < num_eqn; eqn_index++)
        {
            wave_norm_squared += pow(wave[wave_index*num_waves+eqn_index],2);
            left_dot_product += wave[wave_index*num_waves+eqn_index] * wave_left[wave_index*num_waves+eqn_index];
            right_dot_product += wave[wave_index*num_waves+eqn_index] * wave_right[wave_index*num_waves+eqn_index];
        }
        
        if (s[wave_index] < 0.0)
        {
            r = left_dot_product / wave_norm_squared;
        }
        else
        {
            r = right_dot_product / wave_norm_squared;
        }
        
        for (eqn_index=0; eqn_index < num_eqn; eqn_index++)
        {
            wave[wave_index*num_waves+eqn_index] = limiter_func(r) * wave[wave_index*num_waves+eqn_index];
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
