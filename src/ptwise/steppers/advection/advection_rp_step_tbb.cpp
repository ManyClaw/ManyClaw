#include "advection_rp_step_tbb.h"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>

struct advection_rp_step_tbb_body
{
  const real* q;
  const real* aux;
  int nx, ny;
  real* amdq; real* apdq;
  real* wave; real* wave_speeds;

  advection_rp_step_tbb_body(const real* q,  const real* aux,
                             const int nx, const  int ny,
                             real* amdq, real* apdq, real* wave,
                             real* wave_speeds)
    : q(q), aux(aux), nx(nx), ny(ny), amdq(amdq), apdq(apdq),
        wave(wave), wave_speeds(wave_speeds)
  {}

  void operator()(const tbb::blocked_range2d<int>& r) const
  {

    int col, row, idx_left, idx_center, idx_up, idx_out_x, idx_out_y;
    const int num_ghost = advection_rp_grid_params.num_ghost;
    const int num_states = advection_rp_grid_params.num_states;
    const int num_waves = advection_rp_grid_params.num_waves;

    for(row = r.cols().begin(); row < r.cols().end(); ++row) {
      for(col = r.rows().begin(); col < r.rows().end(); ++col) {
        idx_left = col + row*(nx + 2*num_ghost) - 1;
        idx_up = col + (row - 1)*(nx + 2*num_ghost);
        idx_center = idx_left + 1;
        idx_out_x = (col - num_ghost) + (row - num_ghost) * (nx + 1);
        idx_out_y = idx_out_x + ((nx + 1)*(ny + 1));
        advection_rp(q + idx_left*num_states, q + idx_center*num_states,
                     aux, aux, &advection_rp_aux_global,
                     amdq + idx_out_x*num_states, apdq + idx_out_x*num_states,
                     wave + num_waves*num_states*idx_out_x, wave_speeds + num_waves*idx_out_x);
        advection_rp(q + idx_up*num_states, q + idx_center*num_states,
                     aux, aux, &advection_rp_aux_global,
                     amdq + idx_out_y*num_states, apdq + idx_out_y*num_states,
                     wave + num_waves*num_states*idx_out_y, wave_speeds + num_waves*idx_out_y);
      }
    }
  }
};

void advection_rp_step_tbb( const real* q,  const real* aux,
                            const int nx, const  int ny,
                            real* amdq, real* apdq, real* wave,
                            real* wave_speeds)
{
  advection_rp_step_tbb_body body
    (q, aux, nx, ny, amdq, apdq, wave, wave_speeds);

  const int num_ghost = advection_rp_grid_params.num_ghost;
  // note: we use nx+1 and ny+1 here and < in the body (instead of <= in the serial reference)
  tbb::parallel_for(::tbb::blocked_range2d<int,int>(num_ghost, nx+num_ghost+1, num_ghost, ny+num_ghost+1), body);
}

