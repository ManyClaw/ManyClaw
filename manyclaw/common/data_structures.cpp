#include "data_structures.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

Grid::Grid(int nx, int ny) :
  nx(nx), ny(ny)
{
  // TODO make this a function of the encapsulated states
  // Cell centered values
}

State::State(Grid& grid, int num_states, int num_aux, int num_ghost) :
  num_states(num_states), num_aux(num_aux), num_ghost(num_ghost), grid(grid)
{
  const int nx = grid.nx;
  const int ny = grid.ny;
  q.resize((nx+num_ghost*2)*(ny+num_ghost*2)*num_states);
  aux.resize((nx+num_ghost*2)*(ny+num_ghost*2)*num_aux);
}

void Solution::write(int frame, char output_path[])
{
  const char file_prefix[] = "fort.";

  std::ostringstream q_file_path;
  q_file_path << output_path << "/" << file_prefix << "q"
              << std::setfill('0') << std::setw(4) << frame;
  std::cout << "Writing out frame: " << frame << '\n';

  std::fstream q_file;
  q_file.open(q_file_path.str().c_str(), std::ios::out);

  // Write header for fort.q file
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(5) 
              << 1 << "                 grid_number\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(5)
              << 1 << "                 AMR_level\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(5)
              << grid.nx << "                 mx\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(5)
              << grid.ny << "                 my\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(8)
              << std::setw(18) << grid.lower[0] << "    xlow\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(8)
              << std::setw(18) << grid.lower[1] << "    ylow\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(8)
              << std::setw(18) << grid.dx[0] << "    dx\n\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(8)
              << std::setw(18) << grid.dx[1] << "    dy\n\n";

  // Write out q data
  for (int j = state.num_ghost; j < grid.ny + state.num_ghost; j++)
  {
    for(int i = state.num_ghost; i < grid.nx + state.num_ghost; i++)
    {
      for (int m = 0; m < state.num_states; m++)
      {
        q_file << std::setiosflags(std::ios::fixed) 
               << std::setprecision(8)
               << std::setw(16) 
               << state.q[m + i * state.num_states 
                            + j * (2*state.num_ghost + grid.nx)] 
               << " ";
      }
      q_file << "\n";
    }
    q_file << "\n";
  }

  q_file.close();

  // Write out fort.t file
  std::ostringstream t_file_path;
  t_file_path << output_path << "/" << file_prefix << "t"
              << std::setfill('0') << std::setw(4) << frame;

  std::fstream t_file;
  t_file.open(t_file_path.str().c_str(), std::ios::out);

  t_file << std::setiosflags(std::ios::fixed) << std::setprecision(16)
              << std::setw(26) << t << "        time\n";
  t_file << std::setiosflags(std::ios::fixed) << std::setprecision(5) 
              << state.num_states << "                 meqn\n";
  t_file << std::setiosflags(std::ios::fixed) << std::setprecision(5) 
              << 1 << "                 ngrids\n";
  t_file << std::setiosflags(std::ios::fixed) << std::setprecision(5) 
              << 0 << "                 maux\n";
  t_file << std::setiosflags(std::ios::fixed) << std::setprecision(5) 
              << 2 << "                 ndim\n";

  t_file.close();
}

void State::randomize()
{
  randomize_vector(q);
  randomize_vector(aux);
}

Solver::Solver(Solution& solution, int num_waves):
  num_waves(num_waves), solution(solution)
{
  const int nx = solution.grid.nx;
  const int ny = solution.grid.ny;
  const int num_states = solution.state.num_states;
  const int dim = solution.grid.dim;
  // Outputs on interfaces
  amdq.resize((nx+1)*(ny+1)*num_states*dim);
  apdq.resize((nx+1)*(ny+1)*num_states*dim);
  waves.resize((nx+1)*(ny+1)*num_states*num_waves*dim);
  wave_speeds.resize((nx+1)*(ny+1)*4*num_states*dim);
}


void Solver::step(Solution& solution, double dt, 
                  rp_step_t rp_step, updater_t update)
{
  // Note that this all will break if the grid is not uniform!
  real dtdx = dt / solution.grid.dx[0];

  rp_step(&solution.state.q[0], &solution.state.aux[0], 
          solution.grid.num_cells[0], solution.grid.num_cells[1],
          &amdq[0], &apdq[0], &waves[0], &wave_speeds[0]);

  update(&solution.state.q[0], &solution.state.aux[0], 
          solution.grid.num_cells[0], solution.grid.num_cells[1],
          &amdq[0], &apdq[0], &waves[0], &wave_speeds[0],
          num_ghost, solution.state.num_states, dtdx);

  solution.t += dt;
}