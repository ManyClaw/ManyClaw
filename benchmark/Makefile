
CXX=g++
LINK ?= $(CXX)
ifeq ($(CXX), g++)
	OMPCFLAGS = -fopenmp
else ifeq ($(CXX), icc)
	OMPCFLAGS = -openmp
endif
CXXFLAGS=-I ./ -O3 -Wall $(OMPCFLAGS) -I..
LDFLAGS= -ltbb -L../lib -lmanyclaw

DRIVERS = advection \
          euler \
          acoustics_const \
          acoustics_var

all: $(DRIVERS)

advection: advection.o 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

euler: euler.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

acoustics_const: acoustics_const.o 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

acoustics_var: acoustics_var.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

clean:
	rm -rf advection advection.o euler euler.o acoustics_const acoustics_const.o acoustics_var acoustics_var.o
