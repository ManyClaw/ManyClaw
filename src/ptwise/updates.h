
// Upwind updater
void updater_first_order_dimensional_splitting(real* q,
                                               const real* aux,
                                               const int nx,
                                               const int ny,
                                               const real* amdq,
                                               const real* apdq,
                                               const real* wave,
                                               const real* wave_speeds,
                                               const rp_grid_params grid_params
                                               )
{
  int col, row, idx_left, idx_center, idx_up, idx_out_x, idx_out_y;
  const int num_ghost = grid_params.num_ghost;
  const int num_states = grid_params.num_states;
  //  const int num_waves = rp_grid_params.num_waves;

#pragma omp parallel for schedule(runtime) nowait
  for(row = num_ghost; row <= ny + num_ghost; ++row) {
    for(col = num_ghost; col <= nx + num_ghost; ++col) {
      idx_left = col + row*(nx + 2*num_ghost) - 1;
      idx_up = col + (row - 1)*(nx + 2*num_ghost);
      idx_center = idx_left + 1;
      idx_out_x = (col - num_ghost) + (row - num_ghost) * (nx + 1);
      idx_out_y = idx_out_x + ((nx + 1)*(ny + 1));

      for(int state=0; state < num_states; ++state){
        q[idx_left*num_states + state]   -= amdq[idx_out_x*num_states + state];
        q[idx_up*num_states + state]     -= amdq[idx_out_y*num_states + state];
        q[idx_center*num_states + state] -= apdq[idx_out_x*num_states + state];
        q[idx_center*num_states + state] -= apdq[idx_out_y*num_states + state];
      }
    }
  }
}

