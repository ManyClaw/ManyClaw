CXX=g++
CXXFLAGS=-I./ -O3 -Wall
ISPC=ispc
ISPCFLAGS=-O2 --arch=x86-64 --target=sse2

default: driver

.PHONY: dirs clean
.PRECIOUS: objs/simple.h

clean:
	/bin/rm -rf driver *.o

driver: driver.o advection_rp_step_serial.o
	$(CXX) $(CXXFLAGS) -o $@ *.o

advection_rp_step_serial.o: advection_rp_step_serial.h
	$(CXX) $(CXXFLAGS) -o $@ -c advection_rp_step_serial.cpp

#objs/simple.o: simple.cpp objs/simple_ispc.h 
#	$(CXX) $(CXXFLAGS) -c -o $@ $<
#
#objs/%_ispc.h objs/%_ispc.o: %.ispc
#	$(ISPC) $(ISPCFLAGS) $< -o objs/$*_ispc.o -h objs/$*_ispc.h
