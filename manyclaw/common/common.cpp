#include "data_structures.h"
#include "common.h"
#include "timer.h"

void compare_steppers(int nx, int ny, rp_grid_params params,
                      rp_step_t rp_stepper_1, rp_step_t rp_stepper_2)
{
  Grid grid(nx, ny);

  State state(grid, params.num_eqn, params.num_aux, params.num_ghost);
  state.randomize();
  Solution solution(grid, state);

  Solver solver_1(solution, params.num_ghost, params.num_wave);
  Solver solver_2(solution, params.num_ghost, params.num_wave);

  rp_stepper_1(&state.q[0], &state.aux[0], grid.num_cells[0], grid.num_cells[1],
              &solver_1.amdq[0], &solver_1.apdq[0], &solver_1.wave[0],
              &solver_1.wave_speed[0]);
  rp_stepper_2(&state.q[0], &state.aux[0], grid.num_cells[0], grid.num_cells[1],
              &solver_2.amdq[0], &solver_2.apdq[0], &solver_2.wave[0],
              &solver_2.wave_speed[0]);

  std::cout << "  Comparing solution to reference solution...\n";
  std::cout << "    amdq        " << max_error(solver_1.amdq, solver_2.amdq) << "\n";
  std::cout << "    apdq        " << max_error(solver_1.apdq, solver_2.apdq) << "\n";
  std::cout << "    wave       " << max_error(solver_1.wave, solver_2.wave) << "\n";
  std::cout << "    wave_speed " << max_error(solver_1.wave_speed, solver_2.wave_speed) << "\n";
}

double benchmark_stepper(int nx, int ny, rp_grid_params params, rp_step_t rp_stepper)
{

  Grid grid(nx, ny);

  State state(grid, params.num_eqn, params.num_aux, params.num_ghost);
  state.randomize();
  Solution solution(grid, state);
  Solver solver(solution, params.num_ghost, params.num_wave);

  timer t;
  rp_stepper(&state.q[0], &state.aux[0], grid.num_cells[0], grid.num_cells[1],
             &solver.amdq[0], &solver.apdq[0], &solver.wave[0],
             &solver.wave_speed[0]);
  double time = t.elapsed();

  return time;
}


double compare_updates(int nx, int ny, rp_grid_params params, 
                        rp_step_t rp_stepper, updater_t updater)
{
  int index;
  Grid grid(nx, ny);

  State state(grid, params.num_eqn, params.num_aux, params.num_ghost);
  State ref_state(grid, params.num_eqn, params.num_aux, params.num_ghost);

  for (int row = params.num_ghost; row <= ny + params.num_ghost; ++row) {
    for (int col = params.num_ghost; col <= nx + params.num_ghost; ++col) {
      index = col + row * (nx + 2 * params.num_ghost);
      for (int state = 0; state < params.num_eqn; ++state) {
        // q[index * params.num_eqn + state] = 0.0;
        // q[index * params.num_eqn + state] = 0.0;

      }
    }
  }

  state.randomize();
  ref_state.randomize();
 
  Solution solution(grid, state);
  Solver solver(solution, params.num_ghost, params.num_wave);

  rp_stepper(&state.q[0], &state.aux[0], grid.num_cells[0], grid.num_cells[1],
              &solver.amdq[0], &solver.apdq[0], &solver.wave[0],
              &solver.wave_speed[0]);

  updater(&state.q[0], &state.aux[0], grid.num_cells[0], grid.num_cells[1], 
          &solver.amdq[0],  &solver.apdq[0],  &solver.wave[0], 
          &solver.wave_speed[0], params.num_ghost, params.num_eqn, 1.0);

  return max_error(ref_state.q, state.q);
}
