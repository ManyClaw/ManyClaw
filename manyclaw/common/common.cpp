#include "data_structures.h"
#include "common.h"
#include "timer.h"

double benchmark_grid_eval(int nx, int ny, rp_grid_params_t params, 
                           rp_grid_eval_t rp_grid_eval,
                           void* aux_global)
{

  Grid grid(nx, ny);

  State state(grid, params.num_eqn, params.num_aux, params.num_ghost, aux_global);
  state.randomize();
  Solution solution(grid, state);
  Solver solver(solution, params.num_ghost, params.num_wave);

  timer t;
  rp_grid_eval(&state.q[0], &state.aux[0], aux_global, grid.num_cells[0], grid.num_cells[1],
             &solver.amdq[0], &solver.apdq[0], &solver.wave[0],
             &solver.wave_speed[0]);
  double time = t.elapsed();

  return time;
}

