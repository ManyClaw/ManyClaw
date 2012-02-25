include ../variables.inc
include ../rules.inc

all: advection euler


ADVECTION_OBJS = advection_rp_step_serial.o advection_rp_step_serial_cellwise.o advection_rp_step_serial_tiled.o advection_rp_step_omp.o advection_rp_step_tbb.o
 #advection_rp_step_ispc.o advection_rp_step_ispc.h

advection: advection.o
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ $@.o $(LDFLAGS) $(addprefix ../ptwise/steppers/advection/,$(ADVECTION_OBJS)) ../common/common.o ../common/data_structures.o

EULER_OBJS = euler_rp_step_serial.o euler_rp_step_serial_cellwise.o euler_rp_step_serial_tiled.o euler_rp_step_omp.o euler_rp_step_tbb.o
 #euler_rp_step_ispc.o euler_rp_step_ispc.h

euler: euler.o
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ $@.o $(LDFLAGS) $(addprefix ../ptwise/steppers/euler/,$(EULER_OBJS)) ../common/common.o ../common/data_structures.o

clean:
	rm -rf advection advection.o euler euler.o
