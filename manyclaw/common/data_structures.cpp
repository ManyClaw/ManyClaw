#include "data_structures.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

Grid::Grid(int nx, int ny)
{
  // TODO make this a function of the encapsulated states
  // Cell centered values
  num_cells[0] = nx;
  num_cells[1] = ny;
}

State::State(Grid& grid, int num_eqn, int num_aux, int num_ghost, const void* aux_global) :
  num_eqn(num_eqn), num_aux(num_aux), num_ghost(num_ghost), aux_global(aux_global), grid(grid)
{
  const int nx = grid.num_cells[0];
  const int ny = grid.num_cells[1];
  q.resize((nx + num_ghost * 2) * (ny + num_ghost * 2) * num_eqn);
  aux.resize((nx + num_ghost * 2) * (ny + num_ghost * 2) * num_aux);
}

void Solution::write(int frame, std::string output_path)
{
  const char file_prefix[] = "fort.";

  // TODO: Create output directory if not present

  std::ostringstream q_file_path;
  q_file_path << output_path << "/" << file_prefix << "q"
              << std::setfill('0') << std::setw(4) << frame;
  std::cout << "Writing out frame: " << frame << " to " << q_file_path.str() << '\n';

  std::fstream q_file;
  q_file.open(q_file_path.str().c_str(), std::ios::out);

  // Write header for fort.q file
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(5) 
              << 1 << "                 grid_number\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(5)
              << 1 << "                 AMR_level\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(5)
              << grid.num_cells[0] << "                 mx\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(5)
              << grid.num_cells[1] << "                 my\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(8)
              << std::setw(18) << grid.lower[0] << "    xlow\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(8)
              << std::setw(18) << grid.lower[1] << "    ylow\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(8)
              << std::setw(18) << grid.dx[0] << "    dx\n\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(8)
              << std::setw(18) << grid.dx[1] << "    dy\n\n";

  // Write out q data
  FieldIndexer fi(grid.num_cells[0], grid.num_cells[1], state.num_ghost,
                  state.num_eqn);
  for (int row = state.num_ghost; row < grid.num_cells[1] + state.num_ghost; ++row)
  {
    for (int col = state.num_ghost; col < grid.num_cells[0] + state.num_ghost; ++col)
    {
      for (int eqn = 0; eqn < state.num_eqn; ++eqn)
      {
        q_file << std::setiosflags(std::ios::fixed) 
               << std::setprecision(8)
               << std::setw(16) 
               << state.q[eqn + fi.idx(row, col)] 
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
              << state.num_eqn << "                 meqn\n";
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

Solver::Solver()
{
  // Default constructor;
}

Solver::Solver(Solution& solution, int num_ghost, int num_wave):
  num_ghost(num_ghost), num_wave(num_wave)
{
  const int nx = solution.grid.num_cells[0];
  const int ny = solution.grid.num_cells[1];
  const int num_eqn = solution.state.num_eqn;
  //const int dim = solution.grid.dim;
  // Outputs on interfaces
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer wave_efi(nx, ny, num_ghost, num_eqn, num_wave);
  amdq.resize(efi.size());
  apdq.resize(efi.size());
  wave.resize(wave_efi.size());
  wave_speed.resize(wave_efi.size());
}

Solver::Solver(int* num_cells, int num_eqn, int num_ghost, int num_wave):
  num_ghost(num_ghost), num_wave(num_wave)
{
  const int nx = num_cells[0];
  const int ny = num_cells[1];
  //const int dim = solution.grid.dim;
  // Outputs on interfaces
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer wave_efi(nx, ny, num_ghost, num_eqn, num_wave);
  amdq.resize(efi.size());
  apdq.resize(efi.size());
  wave.resize(wave_efi.size());
  wave_speed.resize(wave_efi.size());
}

void Solver::define(int* num_cells, int num_eqn, int num_ghost, int num_wave)
{
  this->num_wave = num_wave;
  this->num_ghost = num_ghost;
  const int nx = num_cells[0];
  const int ny = num_cells[1];

  // Outputs on interfaces
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer wave_efi(nx, ny, num_ghost, num_eqn, num_wave);
  amdq.resize(efi.size());
  apdq.resize(efi.size());
  wave.resize(wave_efi.size());
  wave_speed.resize(wave_efi.size());
}

void Solver::step(Solution& solution, double dt, set_bc_t set_bc, 
                                                 rp_grid_eval_t rp_grid_eval, 
                                                 updater_t update)
{
  // Note that this all will break if the grid is not uniform!
  real dtdx = dt / solution.grid.dx[0];

  set_bc(&solution.state.q[0], &solution.state.aux[0],
        solution.grid.num_cells[0], solution.grid.num_cells[1],
        num_ghost, solution.state.num_eqn);

  rp_grid_eval(&solution.state.q[0], &solution.state.aux[0], 
                solution.state.aux_global, solution.grid.num_cells[0], 
                solution.grid.num_cells[1],
               &amdq[0], &apdq[0], &wave[0], &wave_speed[0]);

  update(&solution.state.q[0], &solution.state.aux[0], 
          solution.grid.num_cells[0], solution.grid.num_cells[1],
          &amdq[0], &apdq[0], &wave[0], &wave_speed[0],
          num_ghost, solution.state.num_eqn, dtdx);

  solution.t += dt;
}
