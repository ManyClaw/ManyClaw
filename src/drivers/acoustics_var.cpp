#include "../many_claw.h"

#include "../ptwise/riemann/acoustics_var_rp.h"
#include "../ptwise/steppers/acoustics_var/acoustics_var_rp_step_serial.h"
#include "../ptwise/steppers/acoustics_var/acoustics_var_rp_step_serial_tiled.h"
#include "../ptwise/steppers/acoustics_var/acoustics_var_rp_step_serial_cellwise.h"
#include "../ptwise/steppers/acoustics_var/acoustics_var_rp_step_tbb.h"
//#include "acoustics_var_rp_step_ispc.h"
#include "../ptwise/steppers/acoustics_var/acoustics_var_rp_step_omp.h"
// TODO add other step headers here

int main(int argc, char ** argv)
{
  int nx = 100, ny = 100;

  if (argc == 3)
  {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
  }

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
    std::cout << "  Benchmark finished in " <<
      1e3 * benchmark_stepper(nx, ny, acoustics_var_rp_grid_params, acoustics_var_rp_steppers[i]) << " ms\n";
  }

  return 0;
}

