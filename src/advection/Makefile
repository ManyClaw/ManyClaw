CXX=g++
CXXFLAGS=-I./ -O3 -Wall
LDFLAGS= #-ltbb
ISPC=ispc
ISPCFLAGS=-O2 --arch=x86-64 --target=sse2 #,sse4-x2,avx-x2
OMPCFLAGS=-fopenmp

default: driver

.PHONY: dirs clean
.PRECIOUS: objs/simple.h

clean:
	/bin/rm -rf driver *.o

driver: driver.o advection_rp_step_serial.o advection_rp_step_omp.o # advection_rp_step_tbb.o advection_rp_step_ispc.o advection_rp_step_ispc.h
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ *.o $(LDFLAGS)

# ISPC must generate the advection_rp_step_ispc.h header before compiling driver
driver.o: driver.cpp #advection_rp_step_ispc.o
	$(CXX) $(CXXFLAGS) -o $@ -c driver.cpp

advection_rp_step_serial.o: advection_rp_step_serial.cpp advection_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c advection_rp_step_serial.cpp

advection_rp_step_tbb.o: advection_rp_step_tbb.cpp advection_rp.h
	$(CXX) $(CXXFLAGS) -o $@ -c advection_rp_step_tbb.cpp

advection_rp_step_ispc.o: advection_rp_step_ispc.ispc advection_rp.h
	$(ISPC) $(ISPCFLAGS) -o $@ -h advection_rp_step_ispc.h advection_rp_step_ispc.ispc

advection_rp_step_omp.o: advection_rp_step_omp.cpp advection_rp.h
	$(CXX) $(CXXFLAGS) $(OMPCFLAGS) -o $@ -c advection_rp_step_omp.cpp

