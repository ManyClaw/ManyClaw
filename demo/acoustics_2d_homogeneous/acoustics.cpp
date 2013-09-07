#include <manyclaw/manyclaw.h>
#include <tbb/task_scheduler_init.h>
#include <string>

#define _USE_MATH_DEFINES

#include <cmath>

int main(int argc, char ** argv)
{
  int nx = 100, ny = 100;

  if (argc == 3)
  {
    nx = std::atoi(argv[1]);
    ny = std::atoi(argv[2]);
  }

  int num_eqn = 3;
  int num_aux = 0;
  int num_ghost = 2;
  int num_wave = 2;
  std::string output_path = "./_output";

  // Initialize solution
  Grid grid(nx, ny);
  grid.lower[0] = -1.0;
  grid.upper[0] =  1.0;
  grid.lower[1] = -1.0;
  grid.upper[1] =  1.0;
  for (int dim = 0; dim < grid.dim; dim++)
    grid.dx[dim] = (grid.upper[dim] - grid.lower[dim]) / grid.num_cells[dim];

  // Physical parameters
  real rho, bulk_modulus, sound_speed, impedence;
  rho = 1.0;
  bulk_modulus = 4.0;
  sound_speed = sqrt(bulk_modulus / rho);
  impedence = rho * sound_speed;
  acoustics_const_rp_aux_global_t aux_global = {rho, bulk_modulus, 
                                                sound_speed, impedence};

  State state(grid, num_eqn, num_aux, num_ghost, &aux_global);
  Solution solution(grid, state);
  solution.t = 0.0;

  // Initialize q  
  const real width = 0.2;
  real x, y, r;
  FieldIndexer fi_q(nx, ny, num_ghost, num_eqn);
  FieldIndexer fi_aux(nx, ny, num_ghost, num_aux);
  for (int row = num_ghost; row < ny + num_ghost; ++row)
  {
    y = grid.lower[1] + (real(row) - 1.5) * grid.dx[1];
    for (int col = num_ghost; col < nx + num_ghost; ++col)
    {
      x = grid.lower[0] + (real(col) - 1.5) * grid.dx[0];
      r = sqrt(pow(x, 2) + pow(y, 2));

      // Background state
      solution.state.q[0 + fi_q.idx(row, col)] = 0.0;
      solution.state.q[1 + fi_q.idx(row, col)] = 0.0;
      solution.state.q[2 + fi_q.idx(row, col)] = 0.0;

      // Perturb within circular radis of width
      if (fabs(r - 0.5) <= width)
        solution.state.q[0 + fi_q.idx(row, col)] = 
                                            1.0 + cos(M_PI * (r - 0.5) / width);
    }
  }
  
  // Write out initial condition
  solution.write(0, output_path);

  // Initialize solver
  Solver solver(solution, num_ghost, num_wave);
  solver.num_ghost = num_ghost;

  // Take multiple steps
  double dt = grid.dx[0] * 0.4 / sound_speed;
  for (int frame = 1; frame <= 20; frame++)
  {
    for (int steps = 0; steps < 10; steps++)
    {
      // Take a single time step
      solver.step(solution, dt, set_zero_order_extrap_BCs, 
                                acoustics_const_rp_grid_eval_void_serial,
                                updater_first_order_dimensional_splitting);
      std::cout << "Solution now at t=" << solution.t << "\n";
    }
    solution.write(frame, output_path);      
  }


  return 0;
}
