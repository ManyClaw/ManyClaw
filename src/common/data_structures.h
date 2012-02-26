#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include "../many_claw.h"
#include "common.h"

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
