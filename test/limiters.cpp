#include <ptwise/limiters.h>

// const int num_eqn, const int num_waves,
//                     const real* wave_left, const real* wave_right,
//                     real (*limiter_func)(real),
//                     real* wave, real* s

#include <iostream>

int main(int argc, char const *argv[])
{
    // Setup test limiter call
    const int num_eqn = 1;
    const int num_waves = 1;

    // State variables
    real wave_right[num_eqn * num_waves];
    real wave_left[num_eqn * num_waves];
    real wave[num_eqn * num_waves];
    real s[num_waves];

    for (int i=0; i < num_eqn; i++)
    {
        for (int m=0; m < num_waves; m++)
        {
            wave_right[i + m * num_waves] = 1.0;
            wave[i + m * num_waves] = 0.5;
            wave_left[i + m * num_waves] = 0.0;
            s[m] = 1.0;
        }
    }

    for (int i=0; i < num_eqn; i++)
        for (int m=0; m < num_waves; m++)
            std::cout << wave[i + m * num_waves] << "\n";

    limiter(1, 1, wave_left, wave_right, beamwarming_limiter, wave, s);

    for (int i=0; i < num_eqn; i++)
        for (int m=0; m < num_waves; m++)
            std::cout << wave[i + m * num_waves] << "\n";

    return 0;
}
