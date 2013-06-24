#include <manyclaw/manyclaw.h>
#include <tbb/task_scheduler_init.h>

int main_kyle(int argc, char ** argv)
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
  // double x, y;
  // for (int j = num_ghost; j < ny + num_ghost; j++)
  // {
  //   y = grid.lower[1] + (double(j) - 1.5) * grid.dx[1];
  //   for (int i = num_ghost; i < nx + num_ghost; i++)
  //   {
  //     x = grid.lower[0] + (double(i) - 1.5) * grid.dx[0];
  //     for (int m = 0; m < num_eqn; m++)
  //     {
  //       if (0.1 < x && x < 0.6 && 0.1 < y && y < 0.6) 
  //         solution.state.q[m + i * num_eqn + j * (2*num_ghost + nx)] = 1.0;
  //       else
  //         solution.state.q[m + i * num_eqn + j * (2*num_ghost + nx)] = 0.1;
  //     }
  //   }
  // }
  for (int j = num_ghost; j < ny + num_ghost; j++)
  {
    // y = grid.lower[1] + (double(j) - 1.5) * grid.dx[1];
    for (int i = num_ghost; i < nx + num_ghost; i++)
    {
      // x = grid.lower[0] + (double(i) - 1.5) * grid.dx[0];
      for (int m = 0; m < num_eqn; m++)
      {
        solution.state.q[m + i * num_eqn + j * (2*num_ghost + nx)] = 0.1;
      }
    }
  }
  solution.state.q[4 * num_eqn + 4 * (2*num_ghost + nx)] = 1.0;
  
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
                                              advection_rp_step_serial_cellwise, 
                                              updater_first_order_dimensional_splitting);
      std::cout << "Solution now at t=" << solution.t << "\n";
    }
    solution.write(frame, "./_output");      
  }


  return 0;
}

int print_solver_q(Solver &solver, int nx, int ny){
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

  return 0;
}

int print_solver(Solver &solver, int nx, int ny){
  print_solver_q(solver, nx, ny);

  std::cout << "solver.waves (size: " << solver.waves.size() << ")\n";
  printf("col  :  ");
  for (int j=0; j < nx+1; ++j)
    printf("%d     ", j);
  printf("\n");
  for (int i=0; i < 2*(ny+1); ++i) {
    printf("row %2d: ", i);
    for (int j=0; j < nx+1; ++j)
      printf("%2.2f  ", solver.waves[j + i*(nx+1)]);
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

  std::cout << "solver.amdq (size: " << solver.amdq.size() << ")\n";
  printf("col  :  ");
  for (int j=0; j < nx+1; ++j)
    printf("%d     ", j);
  printf("\n");
  for (int i=0; i < 2*(ny+1); ++i) {
    printf("row %2d: ", i);
    for (int j=0; j < nx+1; ++j)
      printf("%2.2f  ", solver.amdq[j + i*(nx+1)]);
    printf("\n");
  }
  return 0;
}

int big_test(int argc, char ** argv)
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
  int i= num_ghost, j = num_ghost, m = 0;
  solution.state.q[m + i * num_eqn + j * (2*num_ghost + nx)] = 1.0;

  // Write out initial condition
  solution.write(0, "./_output");

  // Initialize solver
  Solver solver(solution, num_waves, num_ghost);

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

int small_test()
{
  int nx = 5, ny = 5;

  int num_eqn = 1;
  int num_aux = 0;
  int num_ghost = 1;
  int num_waves = 1;

  // Initialize solution
  Grid grid(nx, ny);
  State state(grid, num_eqn, num_aux, num_ghost);
  Solution solution(grid, state);
  Solver solver(solution, num_ghost, num_waves);


  grid.num_cells[0] = nx;
  grid.num_cells[1] = ny;
  grid.lower[0] = 0.0;
  grid.upper[0] = 0.1*nx;
  grid.lower[1] = 0.0;
  grid.upper[1] = 0.1*ny;
  for (int dim = 0; dim < grid.dim; dim++)
    grid.dx[dim] = (grid.upper[dim] - grid.lower[dim]) / grid.num_cells[dim];

  // Initialize q
  int i= num_ghost, j = num_ghost, m = 0;
  solution.state.q[m + i * num_eqn + j * (2*num_ghost + nx)] = 1.0;
  print_solver_q(solver, nx, ny);
  // Take a single time step
  solver.step(solution, grid.dx[0] * 0.9, set_zero_order_extrap_BCs, 
	      advection_rp_step_serial, 
	      updater_first_order_dimensional_splitting);
  print_solver(solver, nx, ny);
  
  return 0;
}

int main(int argc, char ** argv)
{
  //  main_kyle(argc, argv);
  //  big_test(argc, argv);
  small_test();
  return 0;
}
