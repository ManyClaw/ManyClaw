include ../variables.inc

DRIVERS = advection \
          euler \
          acoustics_const \
          acoustics_var

COMMON_OBJS = ../common/common.o ../common/data_structures.o

all: $(DRIVERS)

ADVECTION_STEPS = advection_rp_step_serial.o advection_rp_step_serial_cellwise.o advection_rp_step_serial_tiled.o advection_rp_step_omp.o advection_rp_step_tbb.o
ADVECTION_OBJS = $(addprefix ../ptwise/steppers/advection/,$(ADVECTION_STEPS))

advection: advection.o $(ADVECTION_OBJS) $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ 

EULER_STEPS = euler_rp_step_serial.o euler_rp_step_serial_cellwise.o  euler_rp_step_serial_cellwise_functor.o euler_rp_step_serial_tiled.o euler_rp_step_omp.o euler_rp_step_omp2.o euler_rp_step_tbb.o
EULER_OBJS = $(addprefix ../ptwise/steppers/euler/,$(EULER_STEPS))

euler: euler.o $(EULER_OBJS) $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ 

ACOUSTICS_CONST_STEPS = acoustics_const_rp_step_serial.o acoustics_const_rp_step_serial_cellwise.o acoustics_const_rp_step_serial_tiled.o acoustics_const_rp_step_omp.o acoustics_const_rp_step_tbb.o
ACOUSTICS_CONST_OBJS = $(addprefix ../ptwise/steppers/acoustics_const/,$(ACOUSTICS_CONST_STEPS))

acoustics_const: acoustics_const.o $(ACOUSTICS_CONST_OBJS) $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ 

ACOUSTICS_VAR_STEPS = acoustics_var_rp_step_serial.o acoustics_var_rp_step_serial_cellwise.o acoustics_var_rp_step_serial_tiled.o acoustics_var_rp_step_omp.o acoustics_var_rp_step_tbb.o
ACOUSTICS_VAR_OBJS = $(addprefix ../ptwise/steppers/acoustics_var/,$(ACOUSTICS_VAR_STEPS))

acoustics_var: acoustics_var.o $(ACOUSTICS_VAR_OBJS) $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ 

clean:
	rm -rf advection advection.o euler euler.o acoustics_const acoustics_const.o acoustics_var acoustics_var.o
