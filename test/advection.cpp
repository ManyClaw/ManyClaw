#include <manyclaw/manyclaw.h>
#include <tbb/task_scheduler_init.h>

int main_kyle(int argc, char ** argv)
{
  int nx = 50, ny = 50;

  if (argc == 3)
  {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
  }

  int num_eqn = 1;
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

  State state(grid, num_eqn, num_aux, num_ghost);

  Solution solution(grid, state);
  solution.t = 0.0;

  // Initialize q
  double x, y;
  for (int j = num_ghost; j < ny + num_ghost; j++)
  {
    y = grid.lower[1] + (double(j) - 1.5) * grid.dx[1];
    for (int i = num_ghost; i < nx + num_ghost; i++)
    {
      x = grid.lower[0] + (double(i) - 1.5) * grid.dx[0];
      for (int m = 0; m < num_eqn; m++)
      {
        if (0.1 < x && x < 0.6 && 0.1 < y && y < 0.6) 
          solution.state.q[m + i * num_eqn + j * (2*num_ghost + nx)] = 1.0;
        else
          solution.state.q[m + i * num_eqn + j * (2*num_ghost + nx)] = 0.1;
      }
    }
  }
  
  // Write out initial condition
  solution.write(0, "./_output");

  // Initialize solver
  Solver solver(solution, num_waves);
  solver.num_ghost = num_ghost;

  // Take multiple steps
  for (int frame = 1; frame <= 10; frame++)
  {
    for (int steps = 0; steps < 1; steps++)
    {
      // Take a single time step
      solver.step(solution, grid.dx[0] * 0.9, set_zero_order_extrap_BCs, 
                                              advection_rp_step_serial, 
                                              updater_first_order_dimensional_splitting);
      std::cout << "Solution now at t=" << solution.t << "\n";
    }
    solution.write(frame, "./_output");      
  }


  return 0;
}

int print_solver(Solver &solver, int nx, int ny){
  std::vector<double> &q = solver.solution.state.q;
  std::cout << "solver.q (size: " << q.size() << ")\n";
  printf("col  :  ");
  for (int j=0; j < nx + 2*solver.num_ghost; ++j)
    printf("%d     ", j);
  printf("\n");
  for (int i=0; i < ny + 2*solver.num_ghost; ++i) {
    printf("row %2d: ", i);
    for (int j=0; j < nx + 2*solver.num_ghost; ++j)
      printf("%3.2f  ", q[j + i*(nx+2*solver.num_ghost)]);
    printf("\n");
  }

  std::cout << "solver.apdq (size: " << solver.apdq.size() << ")\n";
  printf("col  :  ");
  for (int j=0; j < nx+1; ++j)
    printf("%d     ", j);
  printf("\n");
  for (int i=0; i < 2*(ny+1); ++i) {
    printf("row %2d: ", i);
    for (int j=0; j < nx+1; ++j)
      printf("%2.2f  ", solver.apdq[j + i*(nx+1)]);
    printf("\n");
  }
  return 0;
}

int main(int argc, char ** argv)
{
  int nx = 10, ny = 10;

  int num_eqn = 1;
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

  State state(grid, num_eqn, num_aux, num_ghost);

  Solution solution(grid, state);
  solution.t = 0.0;

  // Initialize q
  double x, y;
  for (int j = num_ghost; j < ny + num_ghost; j++)
  {
    y = grid.lower[1] + (double(j) - 1.5) * grid.dx[1];
    for (int i = num_ghost; i < nx + num_ghost; i++)
    {
      x = grid.lower[0] + (double(i) - 1.5) * grid.dx[0];
      for (int m = 0; m < num_eqn; m++)
      {
        if (0.1 < x && x < 0.6 && 0.1 < y && y < 0.6) 
          solution.state.q[m + i * num_eqn + j * (2*num_ghost + nx)] = 1.0;
        else
          solution.state.q[m + i * num_eqn + j * (2*num_ghost + nx)] = 0.1;
      }
    }
  }
  
  // Write out initial condition
  solution.write(0, "./_output");

  // Initialize solver
  Solver solver(solution, num_waves);
  solver.num_ghost = num_ghost;

  print_solver(solver, nx, ny);

  // Take multiple steps
  for (int frame = 1; frame <= 3; frame++)
  {
    for (int steps = 0; steps < 1; steps++)
    {
      // Take a single time step
      solver.step(solution, grid.dx[0] * 0.9, set_zero_order_extrap_BCs, 
                                              advection_rp_step_serial, 
                                              updater_first_order_dimensional_splitting);
      std::cout << "Solution now at t=" << solution.t << "\n";
      print_solver(solver, nx, ny);

    }
    solution.write(frame, "./_output");      
  }


  return 0;
}
