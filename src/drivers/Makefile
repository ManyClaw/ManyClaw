include ../variables.inc
include ../rules.inc

all: advection


ADVECTION_OBJS = advection_rp_step_serial.o advection_rp_step_serial_cellwise.o advection_rp_step_serial_tiled.o advection_rp_step_omp.o advection_rp_step_tbb.o
 #advection_rp_step_ispc.o advection_rp_step_ispc.h

advection: advection.o
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ *.o $(LDFLAGS) $(addprefix ../ptwise/steppers/advection/,$(ADVECTION_OBJS)) ../common/common.o ../common/data_structures.o

clean:
	rm -rf advection.o
