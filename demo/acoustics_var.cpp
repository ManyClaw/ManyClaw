#include <manyclaw/manyclaw.h>
#include <tbb/task_scheduler_init.h>

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
  tbb::task_scheduler_init init(num_threads);

  printf("Riemann Solve on %dx%d grid\n", nx, ny);

  const char * acoustics_var_rp_stepper_names[] =
    {
      "serial",
      "serial_tiled",
      "serial_cellwise",
      "TBB",
      "omp"
    };

  rp_step_t acoustics_var_rp_steppers[] =
    {
      acoustics_var_rp_step_serial,
      acoustics_var_rp_step_serial_tiled,
      acoustics_var_rp_step_serial_cellwise,
      acoustics_var_rp_step_tbb,
      acoustics_var_rp_step_omp
      // TODO add other acoustics_var_rp_step functions here
    };

  size_t num_rp_kernels = sizeof(acoustics_var_rp_steppers)/sizeof(rp_step_t);

  for (size_t i = 0; i < num_rp_kernels; i++)
  {
    std::cout << "Testing " << acoustics_var_rp_stepper_names[i] << " Riemann kernel...\n";
    compare_steppers(nx, ny, acoustics_var_rp_grid_params,
                    acoustics_var_rp_step_serial, acoustics_var_rp_steppers[i]);
  }

  return 0;
}

