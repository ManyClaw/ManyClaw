include ../../../variables.inc
include ../../../rules.inc

all: acoustics_const_rp_step_serial.o acoustics_const_rp_step_serial_cellwise.o acoustics_const_rp_step_serial_tiled.o acoustics_const_rp_step_omp.o acoustics_const_rp_step_tbb.o #acoustics_const_rp_step_ispc.o acoustics_const_rp_step_ispc.h

.PHONY: dirs clean
.PRECIOUS: objs/simple.h

clean:
	/bin/rm -rf driver *.o

acoustics_const_driver:
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ *.o $(LDFLAGS)

acoustics_const_rp_step_serial.o: acoustics_const_rp_step_serial.cpp ../../riemann/acoustics_const_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c acoustics_const_rp_step_serial.cpp

acoustics_const_rp_step_serial_tiled.o: acoustics_const_rp_step_serial_tiled.cpp ../../riemann/acoustics_const_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c acoustics_const_rp_step_serial_tiled.cpp

acoustics_const_rp_step_serial_cellwise.o: acoustics_const_rp_step_serial_cellwise.cpp ../../riemann/acoustics_const_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c acoustics_const_rp_step_serial_cellwise.cpp

acoustics_const_rp_step_tbb.o: acoustics_const_rp_step_tbb.cpp ../../riemann/acoustics_const_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c acoustics_const_rp_step_tbb.cpp

acoustics_const_rp_step_ispc.o: acoustics_const_rp_step_ispc.ispc ../../riemann/acoustics_const_rp.h
	$(ISPC) $(ISPCFLAGS) -o $@ -h acoustics_const_rp_step_ispc.h acoustics_const_rp_step_ispc.ispc

acoustics_const_rp_step_omp.o: acoustics_const_rp_step_omp.cpp ../../riemann/acoustics_const_rp.h
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ -c acoustics_const_rp_step_omp.cpp

