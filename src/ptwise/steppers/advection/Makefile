include ../../../variables.inc
include ../../../rules.inc

all: advection_rp_step_serial.o advection_rp_step_serial_cellwise.o advection_rp_step_serial_tiled.o advection_rp_step_omp.o advection_rp_step_tbb.o #advection_rp_step_ispc.o advection_rp_step_ispc.h

.PHONY: dirs clean
.PRECIOUS: objs/simple.h

clean:
	/bin/rm -rf driver *.o

advection_driver:
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ *.o $(LDFLAGS)

advection_rp_step_serial.o: advection_rp_step_serial.cpp ../../riemann/advection_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c advection_rp_step_serial.cpp

advection_rp_step_serial_tiled.o: advection_rp_step_serial_tiled.cpp ../../riemann/advection_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c advection_rp_step_serial_tiled.cpp

advection_rp_step_serial_cellwise.o: advection_rp_step_serial_cellwise.cpp ../../riemann/advection_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c advection_rp_step_serial_cellwise.cpp

advection_rp_step_tbb.o: advection_rp_step_tbb.cpp ../../riemann/advection_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c advection_rp_step_tbb.cpp

advection_rp_step_ispc.o: advection_rp_step_ispc.ispc ../../riemann/advection_rp.h
	$(ISPC) $(ISPCFLAGS) -o $@ -h advection_rp_step_ispc.h advection_rp_step_ispc.ispc

advection_rp_step_omp.o: advection_rp_step_omp.cpp ../../riemann/advection_rp.h
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ -c advection_rp_step_omp.cpp

