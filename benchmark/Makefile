include ../variables.inc

CXXFLAGS += -I../
LDFLAGS =  -L../lib -lmanyclaw $(TBB_LDFLAGS)

DRIVERS = advection \
          advection_var \
          euler \
          acoustics_const \
          acoustics_var

OBJS = benchmark_utils.o

all: $(OBJS) $(DRIVERS) benchit

benchmark_utils.o: benchmark_utils.h

advection: advection.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(OBJS)

advection_var: advection_var.o $(OBJS) 
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(OBJS)

euler: euler.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(OBJS)

acoustics_const: acoustics_const.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(OBJS)

acoustics_var: acoustics_var.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(OBJS)

benchit: $(OBJS) $(DRIVERS)
	@echo "==========================================="
	@echo advection
	./advection
	@echo "==========================================="
	@echo advection_var
	./advection_var
	@echo "==========================================="
	@echo acoustics_const
	./acoustics_const
	@echo "==========================================="
	@echo acoustics_var
	./acoustics_var
	@echo "==========================================="
	@echo euler
	./euler

clean:
	rm -rf advection advection.o advection_var advection_var.o euler euler.o acoustics_const acoustics_const.o acoustics_var acoustics_var.o  $(OBJS)
