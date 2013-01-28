include ../variables.inc
include ../rules.inc

all: advection euler acoustics_const acoustics_var


ADVECTION_OBJS = advection_rp_step_serial.o advection_rp_step_serial_cellwise.o advection_rp_step_serial_tiled.o advection_rp_step_omp.o advection_rp_step_tbb.o
 #advection_rp_step_ispc.o advection_rp_step_ispc.h

advection: advection.o
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ $^ $(LDFLAGS) $(addprefix ../ptwise/steppers/advection/,$(ADVECTION_OBJS)) ../common/common.o ../common/data_structures.o

EULER_OBJS = euler_rp_step_serial.o euler_rp_step_serial_cellwise.o  euler_rp_step_serial_cellwise_functor.o euler_rp_step_serial_tiled.o euler_rp_step_omp.o euler_rp_step_omp2.o euler_rp_step_tbb.o
 #euler_rp_step_ispc.o euler_rp_step_ispc.h

euler: euler.o
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ $^ $(LDFLAGS) $(addprefix ../ptwise/steppers/euler/,$(EULER_OBJS)) ../common/common.o ../common/data_structures.o

ACOUSTICS_CONST_OBJS = acoustics_const_rp_step_serial.o acoustics_const_rp_step_serial_cellwise.o acoustics_const_rp_step_serial_tiled.o acoustics_const_rp_step_omp.o acoustics_const_rp_step_tbb.o
 #acoustics_const_rp_step_ispc.o acoustics_const_rp_step_ispc.h

acoustics_const: acoustics_const.o
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ $^ $(LDFLAGS) $(addprefix ../ptwise/steppers/acoustics_const/,$(ACOUSTICS_CONST_OBJS)) ../common/common.o ../common/data_structures.o

ACOUSTICS_VAR_OBJS = acoustics_var_rp_step_serial.o acoustics_var_rp_step_serial_cellwise.o acoustics_var_rp_step_serial_tiled.o acoustics_var_rp_step_omp.o acoustics_var_rp_step_tbb.o
 #acoustics_var_rp_step_ispc.o acoustics_var_rp_step_ispc.h

acoustics_var: acoustics_var.o
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ $^ $(LDFLAGS) $(addprefix ../ptwise/steppers/acoustics_var/,$(ACOUSTICS_VAR_OBJS)) ../common/common.o ../common/data_structures.o



clean:
	rm -rf advection advection.o euler euler.o acoustics_const acoustics_const.o acoustics_var acoustics_var.o
