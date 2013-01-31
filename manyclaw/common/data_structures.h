#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include <vector>
#include <cmath>
#include <numeric>

typedef double real;

struct rp_grid_params
{
  int num_ghost;
  int num_states;
  int num_aux;
  int num_waves;
};

typedef void (*rp_t)(const real* q_left, const real* q_right,
                  const real* aux_left, const real* aux_right,
                  const void* aux_global,
                     real* amdq, real* apdq, real* wave, real* s);

typedef void (*rp_step_t)(const real* q, const real* aux,
                          const int nx, const int ny, real* amdq,  real* apdq,
                          real* wave, real* wave_speeds);

typedef void (*updater_t)(real* q, const real* aux, const int nx, const int ny,
                                   const real* amdq, const real* apdq,
                                   const real* wave, const real* wave_speeds,
                                   const rp_grid_params grid_params);

template <typename Vector>
void randomize_vector(Vector& v)
{
  // generate numbers in [0,1) deterministically
  real seed = 0.123456789;

  for(size_t i = 0; i < v.size(); i++)
  {
    seed = (314159.26535 * seed);
    seed = seed - floor(seed);
    v[i] = seed;
  }
}

struct Grid
{
  // size info
  int nx;
  int ny;
  static const int dim=2;

  Grid(int nx, int ny);
};

struct State
{
  // size info
  int num_states;
  int num_aux;
  int num_ghost; // not strictly needed but makes life easier for now


  // State and auxilary variables
  std::vector<real> q;
  std::vector<real> aux;

  // Non-owned reference
  Grid &grid;

  State(Grid& grid, int num_states, int num_aux, int num_ghost);
  void randomize();

};

struct Solution
{
  // non-own references
  Grid& grid;
  State& state;

  Solution(Grid& grid, State& state):
    grid(grid), state(state)
  {}
};

struct Solver
{
  // Size info
  int num_ghost;
  int num_waves;

  // Interface data
  std::vector<real> amdq;
  std::vector<real> apdq;
  std::vector<real> waves;
  std::vector<real> wave_speeds;

  // Non-owned references
  Solution& solution;

  Solver(Solution& solution, int num_waves);
};



#endif // DATA_STRUCTURES_H
