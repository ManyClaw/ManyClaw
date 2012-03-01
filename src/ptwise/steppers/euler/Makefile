include ../../../variables.inc
include ../../../rules.inc

all: euler_rp_step_serial.o euler_rp_step_serial_cellwise.o euler_rp_step_serial_tiled.o euler_rp_step_omp.o euler_rp_step_omp2.o euler_rp_step_tbb.o #euler_rp_step_ispc.o euler_rp_step_ispc.h

.PHONY: dirs clean
.PRECIOUS: objs/simple.h

clean:
	/bin/rm -rf driver *.o

euler_driver:
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ *.o $(LDFLAGS)

euler_rp_step_serial.o: euler_rp_step_serial.cpp ../../riemann/euler_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c euler_rp_step_serial.cpp

euler_rp_step_serial_tiled.o: euler_rp_step_serial_tiled.cpp ../../riemann/euler_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c euler_rp_step_serial_tiled.cpp

euler_rp_step_serial_cellwise.o: euler_rp_step_serial_cellwise.cpp ../../riemann/euler_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c euler_rp_step_serial_cellwise.cpp

euler_rp_step_tbb.o: euler_rp_step_tbb.cpp ../../riemann/euler_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c euler_rp_step_tbb.cpp

euler_rp_step_ispc.o: euler_rp_step_ispc.ispc ../../riemann/euler_rp.h
	$(ISPC) $(ISPCFLAGS) -o $@ -h euler_rp_step_ispc.h euler_rp_step_ispc.ispc

euler_rp_step_omp.o: euler_rp_step_omp.cpp ../../riemann/euler_rp.h
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ -c euler_rp_step_omp.cpp

euler_rp_step_omp2.o: euler_rp_step_omp2.cpp ../../riemann/euler_rp.h
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ -c euler_rp_step_omp2.cpp

