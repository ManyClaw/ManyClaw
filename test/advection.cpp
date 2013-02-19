#include <manyclaw/manyclaw.h>
#include <tbb/task_scheduler_init.h>

int main(int argc, char ** argv)
{
  int nx = 10, ny = 10;

  if (argc == 3)
  {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
  }

  int num_eqns = 1;
  int num_aux = 0;
  int num_ghost = 2;
  int num_waves = 1;

  // Initialize solution
  Grid grid(nx, ny);
  grid.num_cells[0] = nx;
  grid.num_cells[1] = ny;
  grid.lower[0] = 0.0;
  grid.upper[0] = 1.0;
  grid.lower[1] = 0.0;
  grid.upper[1] = 1.0;
  for (int dim = 0; dim < grid.dim; dim++)
    grid.dx[dim] = (grid.upper[dim] - grid.lower[dim]) / grid.num_cells[dim];

  State state(grid, num_eqns, num_aux, num_ghost);

  Solution solution(grid, state);
  solution.t = 0.0;

  // for (int j = 0; j < ny + 2*num_ghost; j++)
  //   for (int i = 0; i < nx + 2*num_ghost; i++)
  //     for (int m = 0; m < num_eqns; m++)
  //       solution.state.q[m + i * num_eqns + j * (2*num_ghost + nx)] = 0.0;

  // Initialize q
  double x, y;
  for (int j = num_ghost; j < ny + num_ghost; j++)
  {
    y = grid.lower[1] + (double(j) - 1.5) * grid.dx[1];
    for (int i = num_ghost; i < nx + num_ghost; i++)
    {
      x = grid.lower[0] + (double(i) - 1.5) * grid.dx[0];
      for (int m = 0; m < num_eqns; m++)
      {
        if (0.5 < x && x < 0.75 && 0.5 < y && y < 0.75) 
          solution.state.q[m + i * num_eqns + j * (2*num_ghost + nx)] = 1.0;
        else
          solution.state.q[m + i * num_eqns + j * (2*num_ghost + nx)] = 0.1;
      }
    }
  }
  
  // Write out initial condition
  solution.write(0, "./_output");

  // Initialize solver
  Solver solver(solution, num_waves);
  solver.num_ghost = num_ghost;

  // Take one time step
  solver.step(solution, 0.1);

  solution.write(1, "./_output");

  return 0;
}