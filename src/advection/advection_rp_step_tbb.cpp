#include "advection_rp_step_serial.h"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>

struct advection_rp_step_tbb_body
{
  real* q;
  real* aux;
  int numGhost, numStates, numWaves, nx, ny;
  real* amdq; real* apdq;
  real* wave; real* wave_speeds;

  advection_rp_step_tbb_body
    (real* q,     real* aux,
     int numGhost, int numStates, int numWaves, int nx, int ny, // inputs
     real* amdq,  real* apdq,
     real* wave,  real* wave_speeds)
    : q(q), aux(aux), numGhost(numGhost), numStates(numStates), nx(nx), ny(ny),
      amdq(amdq), apdq(apdq),
      wave(wave), wave_speeds(wave_speeds)
  {}

  void operator()(const tbb::blocked_range2d<int>& r) const
  {
    int x_i, row, left, right;
    for(row = r.cols().begin(); row < r.cols().end(); ++row)
    {
      for(x_i = r.rows().begin(); x_i < r.rows().end(); ++x_i)
      {
        left = x_i + row * nx;
        right = left + 1;
        advection_rp(q + left, q + right, 3,
            aux + left, aux + right,
            amdq + left, apdq + left,
            wave + 2*left, wave_speeds + 2*left);
      }
    }

    int y_i, col;
    for(col = r.rows().begin(); col < r.rows().end(); ++col)
    {
      for(y_i = r.cols().begin(); y_i < r.cols().end(); ++y_i)
      {
        left = y_i * nx + col;
        right = left + nx;
        advection_rp(q + left, q + right, 3,
            aux + left, aux + right,
            amdq + left, apdq + left,
            wave + 2*left, wave_speeds + 2*left);
      }
    }
  }
};

void advection_rp_step_tbb(real* q,     real* aux,
                           int numGhost, int numStates, int numWaves, int nx, int ny, // inputs
                           real* amdq,  real* apdq,
                           real* wave,  real* wave_speeds)
{
  advection_rp_step_tbb_body body
    (q, aux, numGhost, numStates, numWaves, nx, ny,
     amdq, apdq, wave, wave_speeds);

  // note: we use nx+1 and ny+1 here and < in the body (instead of <= in the serial reference)
  tbb::parallel_for(::tbb::blocked_range2d<int,int>(0, nx+1, 0, ny+1), body);
}

