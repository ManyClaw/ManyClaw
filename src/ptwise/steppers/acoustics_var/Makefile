include ../../../variables.inc
include ../../../rules.inc

all: acoustics_var_rp_step_serial.o acoustics_var_rp_step_serial_cellwise.o acoustics_var_rp_step_serial_tiled.o acoustics_var_rp_step_omp.o acoustics_var_rp_step_tbb.o #acoustics_var_rp_step_ispc.o acoustics_var_rp_step_ispc.h

.PHONY: dirs clean
.PRECIOUS: objs/simple.h

clean:
	/bin/rm -rf driver *.o

acoustics_var_driver:
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ *.o $(LDFLAGS)

acoustics_var_rp_step_serial.o: acoustics_var_rp_step_serial.cpp ../../riemann/acoustics_var_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c acoustics_var_rp_step_serial.cpp

acoustics_var_rp_step_serial_tiled.o: acoustics_var_rp_step_serial_tiled.cpp ../../riemann/acoustics_var_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c acoustics_var_rp_step_serial_tiled.cpp

acoustics_var_rp_step_serial_cellwise.o: acoustics_var_rp_step_serial_cellwise.cpp ../../riemann/acoustics_var_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c acoustics_var_rp_step_serial_cellwise.cpp

acoustics_var_rp_step_tbb.o: acoustics_var_rp_step_tbb.cpp ../../riemann/acoustics_var_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c acoustics_var_rp_step_tbb.cpp

acoustics_var_rp_step_ispc.o: acoustics_var_rp_step_ispc.ispc ../../riemann/acoustics_var_rp.h
	$(ISPC) $(ISPCFLAGS) -o $@ -h acoustics_var_rp_step_ispc.h acoustics_var_rp_step_ispc.ispc

acoustics_var_rp_step_omp.o: acoustics_var_rp_step_omp.cpp ../../riemann/acoustics_var_rp.h
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ -c acoustics_var_rp_step_omp.cpp

