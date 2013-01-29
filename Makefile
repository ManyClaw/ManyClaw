include variables.inc

all: lib test benchmark

build:
	$(MAKE) -C src/common
	$(MAKE) -C src/ptwise/steppers/advection
	$(MAKE) -C src/ptwise/steppers/acoustics_const
	$(MAKE) -C src/ptwise/steppers/acoustics_var
	$(MAKE) -C src/ptwise/steppers/euler

lib: build
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o libmanyclaw.so $(shell find src -name "*.o")

test: lib
	$(MAKE) -C test

benchmark: lib
	$(MAKE) -C benchmark

clean:
	$(MAKE) -C src/common clean
	$(MAKE) -C src/ptwise/steppers/advection clean
	$(MAKE) -C src/ptwise/steppers/acoustics_const clean
	$(MAKE) -C src/ptwise/steppers/acoustics_var clean
	$(MAKE) -C src/ptwise/steppers/euler clean
	$(MAKE) -C test clean
	$(MAKE) -C benchmark clean
	rm -rf libmanyclaw.so
