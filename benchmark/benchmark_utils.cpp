#include "benchmark_utils.h"

double benchmark_grid_eval(int nx, int ny, rp_grid_params_t params, 
                           rp_grid_eval_t rp_grid_eval,
                           const void* aux_global)
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

int benchmark(int argc, char **argv, size_t num_kernel, const char *names[], const rp_grid_eval_t evals[],
	      rp_grid_params_t params, const void* aux_global)
{
  int nx = 100, ny = 100;

  if (argc == 3)
  {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
  }

  int num_threads = -1;
  char * descr = getenv("OMP_NUM_THREADS");
  if (descr)
    num_threads = atoi(descr);
  tbb::task_scheduler_init init(num_threads);

  printf("Riemann Solve on %dx%d grid\n", nx, ny);

  for (size_t i = 0; i < num_kernel; i++)
  {
    std::cout << names[i] << "  finished in " <<
      1e3 * benchmark_grid_eval(nx, ny, params, evals[i], aux_global) << " ms\n";
  }

  return 0;
}
