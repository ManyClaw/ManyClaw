#include <manyclaw/manyclaw.h>
#include <tbb/task_scheduler_init.h>
#include <string>

#define _USE_MATH_DEFINES

#include <cmath>

int main(int argc, char ** argv)
{
  int nx = 100, ny = 100, problem = 0;

  if (argc == 3)
  {
    nx = std::atoi(argv[1]);
    ny = std::atoi(argv[2]);
  }
  else if (argc == 4)
  {
    nx = std::atoi(argv[1]);
    ny = std::atoi(argv[2]);
    problem = std::atoi(argv[3]);
  }

  int num_eqn = 4;
  int num_aux = 0;
  int num_ghost = 2;
  int num_wave = 4;
  std::string output_path = "./_output";

  // Physics
  real const gamma = 1.4;
  bool const entropy_fix = false;
  euler_rp_aux_global_t aux_global = {gamma, entropy_fix};

  // Initialize solution
  Grid grid(nx, ny);
  grid.lower[0] = 0.0;
  grid.upper[0] = 1.0;
  grid.lower[1] = 0.0;
  grid.upper[1] = 1.0;
  for (int dim = 0; dim < grid.dim; dim++)
    grid.dx[dim] = (grid.upper[dim] - grid.lower[dim]) / grid.num_cells[dim];

  State state(grid, num_eqn, num_aux, num_ghost, &aux_global);

  Solution solution(grid, state);
  solution.t = 0.0;

  // Initialize q  
  // Quadrant initializations
  real p_init[4], rho_init[4], u_init[4], v_init[4];
  real four_corner_location[2];
  int left_index, right_index;

  // Directionality of Riemann test problems
  if (problem > 0)
  {
    four_corner_location[0] = 0.5;
    four_corner_location[1] = 1.5;

    left_index = 2;
    right_index = 3;
  }
  else
  {
    four_corner_location[0] = 1.5;
    four_corner_location[1] = 0.5;

    left_index = 2;
    right_index = 1;
  }

  if (problem == 0)
  {
    // Quadrants problem
    four_corner_location[0] = 0.8;
    four_corner_location[1] = 0.8;

    // First quadrant (upper right)
    p_init[0] = 1.5;
    rho_init[0] = 1.5;
    u_init[0] = 0.0;
    v_init[0] = 0.0;

    // Second quadrant (upper left)
    p_init[1] = 0.3;
    rho_init[1] = 0.532258064516129;
    u_init[1] = 1.206045378311055;
    v_init[1] = 0.0;

    // Third quadrant (lower left)
    p_init[2] = 0.029032258064516;
    rho_init[2] = 0.137992831541219;
    u_init[2] = 1.206045378311055;
    v_init[2] = 1.206045378311055;

    // Fourth quadrant (lower right)
    p_init[3] = 0.3;
    rho_init[3] = 0.532258064516129;
    u_init[3] = 0.0;
    v_init[3] = 1.206045378311055;
  }
  else if (abs(problem) == 1)
  {
    // Sod test problem
    rho_init[left_index] = 3.0;
    p_init[left_index] = 3.0;
    u_init[left_index] = 0.0;
    v_init[left_index] = 0.0;

    rho_init[right_index] = 1.0;
    p_init[right_index] = 1.0;
    u_init[right_index] = 0.0;
    v_init[right_index] = 0.0;
  }
  else if (abs(problem) == 2)
  {
    // Two shock test problem
    rho_init[left_index] = 1.0;
    u_init[left_index] = 3.0;
    p_init[left_index] = 1.0;
    v_init[left_index] = 0.0;

    rho_init[right_index] = 2.0;
    u_init[right_index] = 1.0;
    p_init[right_index] = 1.0;
    v_init[right_index] = 0.0;
  }

  real x, y;
  int quadrant;
  FieldIndexer fi_q(nx, ny, num_ghost, num_eqn);
  for (int row = num_ghost; row < ny + num_ghost; ++row)
  {
    y = grid.lower[1] + (real(row) - 1.5) * grid.dx[1];
    for (int col = num_ghost; col < nx + num_ghost; ++col)
    {
      x = grid.lower[0] + (real(col) - 1.5) * grid.dx[0];
      
      if (x >= four_corner_location[0] && y >= four_corner_location[1])
        quadrant = 0;
      else if (x < four_corner_location[0] && y >= four_corner_location[1])
        quadrant = 1;
      else if (x < four_corner_location[0] && y < four_corner_location[1])
        quadrant = 2;
      else if (x >= four_corner_location[0] && y < four_corner_location[1])
        quadrant = 3;

      solution.state.q[fi_q.idx(row, col) + 0] = rho_init[quadrant];
      solution.state.q[fi_q.idx(row, col) + 1] = rho_init[quadrant] 
                                                             * u_init[quadrant];
      solution.state.q[fi_q.idx(row, col) + 2] = rho_init[quadrant] 
                                                             * v_init[quadrant];
      solution.state.q[fi_q.idx(row, col) + 3] = p_init[quadrant] / (gamma - 1.0) 
          + 0.5 * rho_init[quadrant] * (  pow(u_init[quadrant],2) 
                                        + pow(v_init[quadrant],2));
    
    }
  }
  
  // Write out initial condition
  solution.write(0, output_path);

  // Initialize solver
  Solver solver(solution, num_ghost, num_wave);
  solver.num_ghost = num_ghost;

  // Take multiple steps
  double dt = grid.dx[0] * 0.04;
  for (int frame = 1; frame <= 20; frame++)
  {
    for (int steps = 0; steps < 100; steps++)
    {
      // Take a single time step
      solver.step(solution, dt, set_zero_order_extrap_BCs, 
                                euler_rp_grid_eval_void_serial,
                                updater_first_order_dimensional_splitting);
      std::cout << "Solution now at t=" << solution.t << "\n";
    }
    solution.write(frame, output_path);      
  }

  return 0;
}
