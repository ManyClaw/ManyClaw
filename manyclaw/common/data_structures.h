#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include <vector>
#include <cmath>
#include <numeric>
#include <string>

typedef double real;

struct rp_grid_params_t
{
  int num_ghost;
  int num_eqn;
  int num_aux;
  int num_wave;
};

typedef void (*set_bc_t)(real* q, real* aux, const int nx, const int ny,
                               const int num_ghost, const int num_eqn);

typedef void (*rp_t)(const real* q_left, const real* q_right,
                     const real* aux_left, const real* aux_right,
                     const void* aux_global, const int direction,
                     real* amdq, real* apdq, real* wave, real* s);

typedef void (*rp_grid_eval_t)(const real* q, const real* aux,
                          const int nx, const int ny, real* amdq,  real* apdq,
                          real* wave, real* wave_speed);

typedef void (*updater_t)(real* q, const real* aux, const int nx, const int ny,
                                   const real* amdq, const real* apdq,
                                   const real* wave, const real* wave_speed,
                                   const int num_ghost, const int num_eqn,
                                   const real dtdx);

template <typename Vector>
void randomize_vector(Vector& v)
{
  // generate numbers in [0,1) deterministically
  real seed = 0.123456789;

  for(unsigned int i = 0; i < v.size(); i++)
  {
    seed = (314159.26535 * seed);
    seed = seed - floor(seed);
    v[i] = seed;
  }
}

struct FieldIndexer
{
  const unsigned nx, ny, num_ghosts, num_eqns;

  FieldIndexer(unsigned nx, unsigned ny, unsigned num_ghosts, unsigned num_eqns)
    : nx(nx), ny(ny), num_ghosts(num_ghosts), num_eqns(num_eqns)
  {}

  inline unsigned idx(int row, int col)
  {return (col + row*(nx + 2*num_ghosts))*num_eqns;}

  inline unsigned up(int row, int col)
  {return (col + (row + 1)*(nx + 2*num_ghosts))*num_eqns;}

  inline unsigned down(int row, int col)
  {return (col + (row - 1)*(nx + 2*num_ghosts))*num_eqns;}

  inline unsigned left(int row, int col)
  {return (col - 1 + row*(nx + 2*num_ghosts))*num_eqns;}

  inline unsigned right(int row, int col)
  {return (col + 1 + row*(nx + 2*num_ghosts))*num_eqns;}

  inline unsigned num_row()
  {return nx + 2 * num_ghosts;}

  inline unsigned num_col()
  {return ny + 2 * num_ghosts;}

  inline unsigned size()
  {return (nx + 2*num_ghosts)*(ny + 2*num_ghosts)*num_eqns;}

  inline unsigned row_size()
  {return (nx + 2*num_ghosts)*num_eqns;}

  inline unsigned col_size()
  {return (ny + 2*num_ghosts)*num_eqns;}
};

// Indexes a field defined on the edges of a grid.
//
// Sample usage:
//   EdgeFieldIndexer(nx, ny, num_ghosts, num_eqns, num_wave);
//   for (int row = 1; row < efi.num_row_edges; ++row) {
//      right = field_data[efi.right_edge(row, col)]
//      left = field_data[efi.left_edge(row, col)]
//   }
//
//  Design Notes:
//  *  num_<foo> is the number of edges based on foo
//  *  <foo>_size is the number of arrays associated with field, so it include
//         striding for multiple equations and wave.
//  *  directions are either transverse (in the direction) or normal (any other direction)
//  *  All members are inline and as const as possible, the hope is we can let 
//         the compiler vectorize any math done on the fields
struct EdgeFieldIndexer
{
  const unsigned nx, ny, num_ghosts, num_eqns, num_wave;

  EdgeFieldIndexer(unsigned nx, unsigned ny, unsigned num_ghosts, 
		   unsigned num_eqns)
    : nx(nx), ny(ny), num_ghosts(num_ghosts), num_eqns(num_eqns), num_wave(1)
  {}

  EdgeFieldIndexer(unsigned nx, unsigned ny, unsigned num_ghosts, 
		   unsigned num_eqns, unsigned num_wave)
    : nx(nx), ny(ny), num_ghosts(num_ghosts), num_eqns(num_eqns), num_wave(num_wave)
  {}

  inline int num_row_edge_normal() const
  {return (nx + 2*num_ghosts - 2);}

  inline int num_row_edge_transverse() const
  {return (nx + 2*num_ghosts - 1);}

  inline int num_row_edge() const
  {return num_row_edge_normal() * num_row_edge_transverse();}

  inline int num_col_edge_normal() const
  {return (ny + 2*num_ghosts - 2);}

  inline int num_col_edge_transverse() const
  {return (ny + 2*num_ghosts - 1);}

  inline int num_col_edge() const
  {return num_col_edge_normal() * num_col_edge_transverse();}

  inline int num_edge() const
  {return num_row_edge() + num_col_edge();}

  inline int row_normal_size() const
  {return num_row_edge_normal()*num_eqns*num_wave;}

  inline int row_transverse_size() const
  {return num_row_edge_transverse()*num_eqns*num_wave;}

  inline int col_normal_size() const
  {return num_col_edge_normal()*num_eqns*num_wave;}

  inline int col_transverse_size() const
  {return num_col_edge_transverse()*num_eqns*num_wave;}

  inline int size() const
  {return num_edge() * num_eqns * num_wave;}

  inline int left_edge(const int row, const int col) const
  {return (col - 1)*num_eqns*num_wave + (row - 1)*col_transverse_size();}

  inline int right_edge(const int row, const int col) const
  {return left_edge(row, col) + 1;}

  inline int down_edge(const int row, const int col) const
  {return (col - 1)*num_eqns*num_wave + (row - 1)*col_normal_size() + size()/2;}

  inline int up_edge(const int row, const int col) const
  {return (col - 1)*num_eqns*num_wave + row*col_normal_size() + size()/2;}
};

struct Grid
{
  // size info
  static const int dim=2;
  int num_cells[dim];
  double dx[dim];
  double lower[dim];
  double upper[dim];

  Grid(int nx, int ny);
};

struct State
{
  // size info
  int num_eqn;
  int num_aux;
  int num_ghost; // not strictly needed but makes life easier for now


  // State and auxilary variables
  std::vector<real> q;
  std::vector<real> aux;

  // Non-owned reference
  Grid &grid;

  State(Grid& grid, int num_eqn, int num_aux, int num_ghost);
  void randomize();

};

struct Solution
{
  // non-own references
  Grid& grid;
  State& state;
  double t;

  Solution(Grid& grid, State& state):
    grid(grid), state(state), t(0.0)
  {}

  void write(int frame, std::string output_path);
};

struct Solver
{
  // Size info
  int num_ghost;
  int num_wave;

  // Interface data
  std::vector<real> amdq;
  std::vector<real> apdq;
  std::vector<real> wave;
  std::vector<real> wave_speed;

  // Non-owned references
  //Solution& solution;

  //TODO: Deprecate this constructor
  Solver(Solution& solution, int num_ghost, int num_wave);

  Solver(int* num_cells, int num_eqn, int num_ghost, int num_wave);
  
  void step(Solution& solution, double dt, set_bc_t set_bc, rp_grid_eval_t rp_grid_eval, updater_t update);

  inline real cfl(int nx, int ny, int mbc, int meqn, int mwaves, real dtdx)
  {
    EdgeFieldIndexer efi(nx, ny, mbc, meqn, mwaves);
    real cfl = 0.0;
    for (int row=mbc; row < ny + mbc + 1; ++row)
    {
      for (int col=mbc; col < nx + mbc + 1; ++col)
      {
        for (int wave=0; wave < num_wave; ++wave)
        {
          cfl = fabs(wave_speed[efi.left_edge(row, col) + wave]) * dtdx;
          cfl = fmax(cfl, fabs(wave_speed[efi.down_edge(row, col) + wave]) * dtdx);
        }
      }
    }
    return cfl;
  }
};



#endif // DATA_STRUCTURES_H
