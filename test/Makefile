include ../variables.inc

CXXFLAGS += -I ..
LDFLAGS += -L../ -lmanyclaw

DRIVERS = advection \
          euler \
          acoustics_const \
          acoustics_var \
          limiters

COMMON_OBJS = ../src/common/common.o ../src/common/data_structures.o

all: $(DRIVERS)

advection: advection.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< 

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
