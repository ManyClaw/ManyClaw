#include "advection_rp_step_serial.h"

#include <iostream>
#include <tbb/blocked_range2d.h>

//void advection_rp_step_serial(real* q,     real* aux,
//                              int numGhost, int numStates, int numWaves, int nx, int ny, // inputs
//                              real* amdq,  real* apdq,
//                              real* wave,  real* wave_speeds)
//{
//  ::tbb::parallel_for(::tbb::blocked_range<Size>(0,n), for_each_detail::make_body<Size>(first,f));
//
//  int x_i, row, left, right;
//  for(row = 0; row <= ny; ++row)
//  {
//    for(x_i = 0; x_i <= nx; ++x_i)
//    {
//      left = x_i + row * nx;
//      right = left + 1;
//      advection_rp(q + left, q + right, 3,
//                   aux + left, aux + right,
//                   amdq + left, apdq + left,
//                   wave + 2*left, wave_speeds + 2*left);
//    }
//  }
//
//  int y_i, col;
//  for(col = 0; col <= nx; ++col)
//  {
//    for(y_i = 0; y_i <= ny; ++y_i)
//    {
//      left = y_i * nx + col;
//      right = left + nx;
//      advection_rp(q + left, q + right, 3,
//                   aux + left, aux + right,
//                   amdq + left, apdq + left,
//                   wave + 2*left, wave_speeds + 2*left);
//    }
//  }
//}

void advection_rp_step_tbb(real* q,     real* aux,
                           int numGhost, int numStates, int numWaves, int nx, int ny, // inputs
                           real* amdq,  real* apdq,
                           real* wave,  real* wave_speeds)
{
  std::cout << "hello from tbb" << std::endl;
}

