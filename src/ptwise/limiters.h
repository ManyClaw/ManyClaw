#ifndef LIMITERS_H
#define LIMITERS_H

#include <math.h>

typedef double real;

// Limit wave argument based on wave_left, wave_right, and s
// This is the monotone centered limiter
inline void limiter(const int num_eqn, const int num_waves,
                    const real* wave_left, const real* wave_right,
                    void *limiter_func,
                    real* wave, real* s)
{
    // Work space
    int eqn_index, wave_index;
    real wave_norm_squared, left_dot_product, right_dot_product;


    for (wave_index = 0; wave_index < num_waves; wave_index++)
    {
        for (eqn_index=0; eqn_index < num_eqn; eqn_index++)
        {
            wave_norm_squared += wave[wave_index,eqn_index]**2;
            left_dot_product += wave[wave_index, eqn_index] * wave_left[wave_index, eqn_index];
            right_dot_product += wave[wave_index, eqn_index] * wave_right[wave_index, eqn_index];
        }

        if (s[wave_index] < 0.0)
        {
            for ()
        }
    }
}

inline void mc_limiter
