#include <manyclaw/manyclaw.h>
#include <tbb/task_scheduler_init.h>

// TODO add other step headers here

int main(int argc, char ** argv)
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
  tbb::task_scheduler_init tsi;
  tsi.initialize(num_threads);

  printf("Riemann Solve on %dx%d grid\n", nx, ny);


  for (size_t i = 0; i < num_advection_rp_grid_eval_kernels; i++)
  {
    std::cout << advection_rp_grid_eval_names[i] << "  finished in " <<
      1e3 * benchmark_grid_eval(nx, ny, advection_rp_grid_params, advection_rp_grid_evals[i]) << " ms\n";
  }
  return 0;
}

