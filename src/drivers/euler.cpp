#include "../many_claw.h"
#include <tbb/task_scheduler_init.h>


#include "../ptwise/riemann/euler_rp.h"
#include "../ptwise/steppers/euler/euler_rp_step_serial.h"
#include "../ptwise/steppers/euler/euler_rp_step_serial_tiled.h"
#include "../ptwise/steppers/euler/euler_rp_step_serial_cellwise.h"
#include "../ptwise/steppers/euler/euler_rp_step_tbb.h"
//#include "euler_rp_step_ispc.h"
#include "../ptwise/steppers/euler/euler_rp_step_omp.h"
#include "../ptwise/steppers/euler/euler_rp_step_omp2.h"
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
  tbb::task_scheduler_init init(num_threads);

  printf("Riemann Solve on %dx%d grid\n", nx, ny);

  const char * euler_rp_stepper_names[] =
    {
      "serial",
      "serial_tiled",
      "serial_cellwise",
      "TBB",
      "omp",
      "omp2"
    };

  rp_step_t euler_rp_steppers[] =
    {
      euler_rp_step_serial,
      euler_rp_step_serial_tiled,
      euler_rp_step_serial_cellwise,
      euler_rp_step_tbb,
      euler_rp_step_omp,
      euler_rp_step_omp2
      // TODO add other euler_rp_step functions here
    };

  size_t num_rp_kernels = sizeof(euler_rp_steppers)/sizeof(rp_step_t);

  for (size_t i = 0; i < num_rp_kernels; i++)
  {
    std::cout << "Testing " << euler_rp_stepper_names[i] << " Riemann kernel...\n";
    compare_steppers(nx, ny, euler_rp_grid_params,
                    euler_rp_step_serial, euler_rp_steppers[i]);
    std::cout << "  Benchmark finished in " <<
      1e3 * benchmark_stepper(nx, ny, euler_rp_grid_params, euler_rp_steppers[i]) << " ms\n";
  }

  return 0;
}

