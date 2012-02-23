CXX=g++
CXXFLAGS=-I./ -O3 -Wall -I/Users/aterrel/workspace/tbb30_20110704oss/include
LDFLAGS= -L/Users/aterrel/workspace/tbb30_20110704oss/lib -ltbb
ISPC=ispc
ISPCFLAGS=-O2 --arch=x86-64 --target=sse2 #,sse4-x2,avx-x2
OMPCFLAGS=-fopenmp

default: advection_driver

.PHONY: dirs clean
.PRECIOUS: objs/simple.h

clean:
	/bin/rm -rf driver *.o

advection_driver: advection_driver.o advection_rp_step_serial.o #advection_rp_step_serial_tiled.o advection_rp_step_serial_cellwise.o advection_rp_step_omp.o advection_rp_step_tbb.o #advection_rp_step_ispc.o advection_rp_step_ispc.h
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ *.o $(LDFLAGS)

# ISPC must generate the advection_rp_step_ispc.h header before compiling driver
advection_driver.o: advection_driver.cpp #advection_rp_step_ispc.o
	$(CXX) $(CXXFLAGS) -o $@ -c advection_driver.cpp

advection_rp_step_serial.o: advection_rp_step_serial.cpp ../ptwise/riemann/advection_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c advection_rp_step_serial.cpp

advection_rp_step_serial_tiled.o: advection_rp_step_serial_tiled.cpp ../ptwise/riemann/advection_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c advection_rp_step_serial_tiled.cpp

advection_rp_step_serial_cellwise.o: advection_rp_step_serial_cellwise.cpp ../ptwise/riemann/advection_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c advection_rp_step_serial_cellwise.cpp

advection_rp_step_tbb.o: advection_rp_step_tbb.cpp ../ptwise/riemann/advection_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c advection_rp_step_tbb.cpp

advection_rp_step_ispc.o: advection_rp_step_ispc.ispc ../ptwise/riemann/advection_rp.h
	$(ISPC) $(ISPCFLAGS) -o $@ -h advection_rp_step_ispc.h advection_rp_step_ispc.ispc

advection_rp_step_omp.o: advection_rp_step_omp.cpp ../ptwise/riemann/advection_rp.h
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ -c advection_rp_step_omp.cpp

