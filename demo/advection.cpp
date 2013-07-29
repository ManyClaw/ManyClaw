#include <manyclaw/manyclaw.h>
#include <tbb/task_scheduler_init.h>
#include <string>

int main(int argc, char ** argv)
{
  int nx = 10, ny = 10;

  if (argc == 3)
  {
    nx = std::atoi(argv[1]);
    ny = std::atoi(argv[2]);
  }

  int num_eqn = 1;
  int num_aux = 0;
  int num_ghost = 2;
  int num_wave = 1;
  std::string output_path = "./_output";

  // Initialize solution
  Grid grid(nx, ny);
  grid.lower[0] = 0.0;
  grid.upper[0] = 1.0;
  grid.lower[1] = 0.0;
  grid.upper[1] = 1.0;
  for (int dim = 0; dim < grid.dim; dim++)
    grid.dx[dim] = (grid.upper[dim] - grid.lower[dim]) / grid.num_cells[dim];

  advection_rp_aux_global_t aux_global = {{1,0.5}};
  State state(grid, num_eqn, num_aux, num_ghost, &aux_global);

  Solution solution(grid, state);
  solution.t = 0.0;

  // Initialize q  
  real x, y;
  FieldIndexer fi_q(nx, ny, num_ghost, num_eqn);
  FieldIndexer fi_aux(nx, ny, num_ghost, num_aux);
  for (int row = num_ghost; row < ny + num_ghost; ++row)
  {
    y = grid.lower[1] + (real(row) - 1.5) * grid.dx[1];
    for (int col = num_ghost; col < nx + num_ghost; ++col)
    {
      x = grid.lower[0] + (real(col) - 1.5) * grid.dx[0];
      if (0.1 < x && x < 0.6 && 0.1 < y && y < 0.6) 
        solution.state.q[0 + fi_q.idx(row, col)] = 1.0;
      else
        solution.state.q[0 + fi_q.idx(row, col)] = 0.1;
    }
  }
  
  // Write out initial condition
  solution.write(0, output_path);

  // Initialize solver
  Solver solver(solution, num_ghost, num_wave);
  solver.num_ghost = num_ghost;

  // Take multiple steps
  for (int frame = 1; frame <= 10; frame++)
  {
    for (int steps = 0; steps < 10; steps++)
    {
      // Take a single time step
      solver.step(solution, grid.dx[0] * 0.4, set_zero_order_extrap_BCs, 
                                              advection_rp_grid_eval_void,
                                              updater_first_order_dimensional_splitting);
      std::cout << "Solution now at t=" << solution.t << "\n";
    }
    solution.write(frame, output_path);      
  }


  return 0;
}
