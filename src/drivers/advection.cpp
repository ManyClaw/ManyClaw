#include "../many_claw.h"
#include <tbb/task_scheduler_init.h>

#include "../ptwise/riemann/advection_rp.h"
#include "../ptwise/steppers/advection/advection_rp_step_serial.h"
#include "../ptwise/steppers/advection/advection_rp_step_serial_tiled.h"
#include "../ptwise/steppers/advection/advection_rp_step_serial_cellwise.h"
#include "../ptwise/steppers/advection/advection_rp_step_tbb.h"
//#include "advection_rp_step_ispc.h"
#include "../ptwise/steppers/advection/advection_rp_step_omp.h"
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

  const char * advection_rp_stepper_names[] =
    {
      "serial",
      "serial_tiled",
      "serial_cellwise",
      "TBB",
      "omp"
    };

  rp_step_t advection_rp_steppers[] =
    {
      advection_rp_step_serial,
      advection_rp_step_serial_tiled,
      advection_rp_step_serial_cellwise,
      advection_rp_step_tbb,
      advection_rp_step_omp
      // TODO add other advection_rp_step functions here
    };

  size_t num_rp_kernels = sizeof(advection_rp_steppers)/sizeof(rp_step_t);

  for (size_t i = 0; i < num_rp_kernels; i++)
  {
    std::cout << "Testing " << advection_rp_stepper_names[i] << " Riemann kernel...\n";
    compare_steppers(nx, ny, advection_rp_grid_params,
                    advection_rp_step_serial, advection_rp_steppers[i]);
    std::cout << "  Benchmark finished in " <<
      1e3 * benchmark_stepper(nx, ny, advection_rp_grid_params, advection_rp_steppers[i]) << " ms\n";
  }

  return 0;
}

