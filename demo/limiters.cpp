#include <manyclaw/ptwise/limiters.h>

// const int num_eqn, const int num_waves,
//                     const real* wave_left, const real* wave_right,
//                     real (*limiter_func)(real),
//                     real* wave, real* s

#include <iostream>
#include <cmath>

double random_number(double seed)
{
    seed = (314159.26535 * seed);
    return seed - floor(seed);
}

int main(int argc, char const *argv[])
{
    // Setup test limiter call
    const int num_eqn = 1;
    const int num_waves = 1;
    const int num_cells = 100;
    const int num_ghost = 2;

    double x;
    double dx = 2.0 / double(num_cells);
    double seed = 0.12345567;

    // State variables
    double *wave = new double[num_eqn * num_waves * (num_cells + 2 * num_ghost)];
    double *s = new double[num_waves * (num_cells + 2 * num_ghost)];

    for (int i=0; i < num_cells + 2 * num_ghost; i++)
    {
        x = -1.0 + (i - 1.5) * dx;
        for (int m=0; m < num_waves; m++)
        {
            for (int eqn_index = 0; eqn_index < num_eqn; eqn_index++)
            {
                // wave[eqn_index + m * num_eqn + i * num_eqn * num_waves] = 
                //         1.0 - pow(x, 2);
                if (x < 0.0)
                    wave[eqn_index + m * num_eqn + i * num_eqn * num_waves] = 0.0;
                else
                    wave[eqn_index + m * num_eqn + i * num_eqn * num_waves] = 1.0;
                        
            }
            seed = random_number(seed);
            s[m + i * num_waves] = seed;
        }
    }

    std::cout << "[";
    for (int i=0; i < num_cells; i++)
        for (int m=0; m < num_waves; m++)
            for (int eqn_index = 0; eqn_index < num_eqn; eqn_index++)
                std::cout << wave[eqn_index + m * num_eqn + i * num_eqn * num_waves] << ", ";
            std::cout << "]\n";
        std::cout << "==================\n";

    limit_waves(num_eqn, num_waves, num_ghost, num_cells, wave, s, beamwarming_limiter);

    std::cout << "[";
    for (int i=0; i < num_cells; i++)
        for (int m=0; m < num_waves; m++)
            for (int eqn_index = 0; eqn_index < num_eqn; eqn_index++)
                std::cout << wave[eqn_index + m * num_eqn + i * num_eqn * num_waves] << ", ";
            std::cout << "]\n";
        std::cout << "==================\n";

    return 0;
}
