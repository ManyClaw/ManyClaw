#include "common.h"

void compare_steppers(int nx, int ny, rp_grid_params params,
                      rp_step_t rp_stepper_1, rp_step_t rp_stepper_2)
{
  Grid grid(nx, ny);

  State state(grid, params.num_states, params.num_aux, params.num_ghost);
  state.randomize();
  Solution solution(grid, state);

  Solver solver_1(solution, params.num_waves);
  Solver solver_2(solution, params.num_waves);

  rp_stepper_1(&state.q[0], &state.aux[0], grid.nx, grid.ny,
              &solver_1.amdq[0], &solver_1.apdq[0], &solver_1.waves[0],
              &solver_1.wave_speeds[0]);
  rp_stepper_2(&state.q[0], &state.aux[0], grid.nx, grid.ny,
              &solver_2.amdq[0], &solver_2.apdq[0], &solver_2.waves[0],
              &solver_2.wave_speeds[0]);

  std::cout << "  Comparing solution to reference solution...\n";
  std::cout << "    amdq        " << max_error(solver_1.amdq, solver_2.amdq) << "\n";
  std::cout << "    apdq        " << max_error(solver_1.apdq, solver_2.apdq) << "\n";
  std::cout << "    waves       " << max_error(solver_1.waves, solver_2.waves) << "\n";
  std::cout << "    wave_speeds " << max_error(solver_1.wave_speeds, solver_2.wave_speeds) << "\n";
}

double benchmark_stepper(int nx, int ny, rp_grid_params params, rp_step_t rp_stepper)
{

  Grid grid(nx, ny);

  State state(grid, params.num_states, params.num_aux, params.num_ghost);
  state.randomize();
  Solution solution(grid, state);
  Solver solver(solution, params.num_waves);

  timer t;
  rp_stepper(&state.q[0], &state.aux[0], grid.nx, grid.ny,
             &solver.amdq[0], &solver.apdq[0], &solver.waves[0],
             &solver.wave_speeds[0]);
  double time = t.elapsed();

  return time;
}
