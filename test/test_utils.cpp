#include "test_utils.h"

void compare_arrays(const real *expected, const real *actual,
                    int size, bool& result, std::ostringstream& msg)
{
  std::vector<int> bad_loc;
  for (int i=0; i < size; ++i){
    if (std::fabs(expected[i] - actual[i]) > 1e-8){
      bad_loc.push_back(i);
    }
  }

  if (bad_loc.size() > 0) {
    msg << std::endl;
    for(size_t i=0; i < bad_loc.size(); ++i){
      int idx = bad_loc[i];
      msg << std::scientific
           << "  array[" << idx
           << "] (" << actual[idx] << ") != expected[" << idx
           << "] (" << expected[idx] << ")\n";
    }
    msg << "  Num_wrong: " << bad_loc.size() << " of " << size;
    result = false;
  } else {
    result = true;
  }
}

::testing::AssertionResult ArraysMatch(const real *actual, const real *expected, int size){
  bool result;
  std::ostringstream msg;
  compare_arrays(expected, actual, size, result, msg);
  if (!result){
    return ::testing::AssertionFailure() << msg.str();
  }
  return ::testing::AssertionSuccess();
}

::testing::AssertionResult GridEvalsMatch(int nx, int ny, rp_grid_params_t params,
                                          rp_grid_eval_t rp_grid_eval_1, const char* rp_grid_eval_name_1, 
                                          rp_grid_eval_t rp_grid_eval_2, const char* rp_grid_eval_name_2, 
                                          const void* aux_global)
{
  bool results[4];
  std::ostringstream msg;
  Grid grid(nx, ny);

  State state(grid, params.num_eqn, params.num_aux, params.num_ghost, aux_global);
  state.randomize();
  Solution solution(grid, state);

  Solver solver_1(solution, params.num_ghost, params.num_wave);
  Solver solver_2(solution, params.num_ghost, params.num_wave);

  rp_grid_eval_1(&state.q[0], &state.aux[0], state.aux_global,
                 grid.num_cells[0], grid.num_cells[1],
                 &solver_1.amdq[0], &solver_1.apdq[0], &solver_1.wave[0],
                 &solver_1.wave_speed[0]);
  rp_grid_eval_2(&state.q[0], &state.aux[0], state.aux_global, 
                 grid.num_cells[0], grid.num_cells[1],
                 &solver_2.amdq[0], &solver_2.apdq[0], &solver_2.wave[0],
                 &solver_2.wave_speed[0]);

  compare_arrays(&solver_1.amdq[0], &solver_2.amdq[0], solver_1.amdq.size(), results[0], msg);
  compare_arrays(&solver_1.apdq[0], &solver_2.apdq[0], solver_1.apdq.size(), results[1], msg);
  compare_arrays(&solver_1.wave[0], &solver_2.wave[0], solver_1.wave.size(), results[2], msg);
  compare_arrays(&solver_1.wave_speed[0], &solver_2.wave_speed[0], solver_1.wave_speed.size(),
                 results[3], msg);
  if (!(results[0] && results[1] && results[2] && results[3]))
    return ::testing::AssertionFailure() << "Error comparing " << rp_grid_eval_name_1 << " and " 
					 << rp_grid_eval_name_2 << "\n" << msg.str();
  return ::testing::AssertionSuccess();
}


double compare_updates(int nx, int ny, rp_grid_params_t params, 
                       rp_grid_eval_t rp_grid_eval, updater_t updater, void* aux_global)
{
  //  int index;
  Grid grid(nx, ny);

  State state(grid, params.num_eqn, params.num_aux, params.num_ghost, aux_global);
  State ref_state(grid, params.num_eqn, params.num_aux, params.num_ghost, aux_global);

  for (int row = params.num_ghost; row <= ny + params.num_ghost; ++row) {
    for (int col = params.num_ghost; col <= nx + params.num_ghost; ++col) {
      //  index = col + row * (nx + 2 * params.num_ghost);
      for (int state = 0; state < params.num_eqn; ++state) {
        // q[index * params.num_eqn + state] = 0.0;
        // q[index * params.num_eqn + state] = 0.0;

      }
    }
  }

  state.randomize();
  ref_state.randomize();
 
  Solution solution(grid, state);
  Solver solver(solution, params.num_ghost, params.num_wave);

  rp_grid_eval(&state.q[0], &state.aux[0], state.aux_global, grid.num_cells[0], grid.num_cells[1],
              &solver.amdq[0], &solver.apdq[0], &solver.wave[0],
              &solver.wave_speed[0]);

  updater(&state.q[0], &state.aux[0], grid.num_cells[0], grid.num_cells[1], 
          &solver.amdq[0],  &solver.apdq[0],  &solver.wave[0], 
          &solver.wave_speed[0], params.num_ghost, params.num_eqn, 1.0);

  return max_error(ref_state.q, state.q);
}
