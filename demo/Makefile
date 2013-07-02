include ../variables.inc

CXXFLAGS += -I..
LDFLAGS = $(TBB_LDFLAGS) -L../lib -lmanyclaw 

DRIVERS = advection
          # euler \
          # acoustics_const \
          # acoustics_var \
          # limiters

all: $(DRIVERS)

advection: advection.o
	$(LINK) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

euler: euler.o 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

acoustics_const: acoustics_const.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

acoustics_var: acoustics_var.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

limiters: limiters.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

clean:
	rm -rf advection advection.o euler euler.o acoustics_const acoustics_const.o acoustics_var acoustics_var.o limiters limiters.o

clobber: clean
	$(MAKE) -C ../ clean
