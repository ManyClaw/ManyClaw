include ../variables.inc

CXXFLAGS += -I../
LDFLAGS =  -L../lib -lmanyclaw $(TBB_LDFLAGS)

DRIVERS = advection 
#\
#          euler \
#          acoustics_const \
#          acoustics_var

all: $(DRIVERS)

advection: advection.o 
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

euler: euler.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) 

acoustics_const: acoustics_const.o 
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

acoustics_var: acoustics_var.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -rf advection advection.o euler euler.o acoustics_const acoustics_const.o acoustics_var acoustics_var.o
