#include "advection_rp_step_serial.h"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>

struct advection_rp_step_tbb_body
{
  real* q;
  real* aux;
  int num_ghost, num_states, num_waves, nx, ny;
  real* amdq; real* apdq;
  real* wave; real* wave_speeds;

  advection_rp_step_tbb_body
    (real* q,     real* aux,
     int num_ghost, int num_states, int num_waves, int nx, int ny, // inputs
     real* amdq,  real* apdq,
     real* wave,  real* wave_speeds)
      : q(q), aux(aux), num_ghost(num_ghost), num_states(num_states),
        num_waves(num_waves), nx(nx), ny(ny), amdq(amdq), apdq(apdq),
        wave(wave), wave_speeds(wave_speeds)
  {}

  void operator()(const tbb::blocked_range2d<int>& r) const
  {
    int col, row, idx_left, idx_right, idx_up, idx_down, idx_out;
    const int num_aux = 2;

    for(row = r.cols().begin(); row < r.cols().end(); ++row)
    {
      for(col = r.rows().begin(); col < r.rows().end(); ++col)
      {
        idx_left = col + row*(nx + 2*num_ghost) - 1;
        idx_right = idx_left + 1;
        idx_out = (col - num_ghost) + (row - num_ghost) * (nx + 1);
        advection_rp(q + idx_left*num_states, q + idx_right*num_states, num_states,
                     aux + idx_left*num_aux, aux + idx_right*num_aux,
                     amdq + idx_out*num_states, apdq + idx_out*num_states,
                     wave + num_waves*num_states*idx_out, wave_speeds + num_waves*idx_out);
      }
    }

    for(col = r.rows().begin(); col < r.rows().end(); ++col)
    {
      for(row = r.cols().begin(); row < r.cols().end(); ++row)
      {
      idx_up = col + (row - 1)*(nx + 2*num_ghost);
      idx_down = nx + 2*num_ghost;
      idx_out = (col - num_ghost) + (row - num_ghost) * (nx + 1) + ((nx + 1)*(ny + 1));
      advection_rp(q + idx_up*num_states, q + idx_down*num_states, num_states,
                   aux + idx_up*num_aux, aux + idx_down*num_aux,
                   amdq + idx_out*num_states, apdq + idx_out*num_states,
                   wave + num_waves*num_states*idx_out, wave_speeds + num_waves*idx_out);
      }
    }
  }
};

void advection_rp_step_tbb(real* q,     real* aux,
                           int num_ghost, int num_states, int num_waves, int nx, int ny, // inputs
                           real* amdq,  real* apdq,
                           real* wave,  real* wave_speeds)
{
  advection_rp_step_tbb_body body
    (q, aux, num_ghost, num_states, num_waves, nx, ny,
     amdq, apdq, wave, wave_speeds);

  // note: we use nx+1 and ny+1 here and < in the body (instead of <= in the serial reference)
  tbb::parallel_for(::tbb::blocked_range2d<int,int>(num_ghost, nx+num_ghost+1, num_ghost, ny+num_ghost+1), body);
}

