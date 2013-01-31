#include "data_structures.h"

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

