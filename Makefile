include variables.inc

all: lib test benchmark

build:
	$(MAKE) -C manyclaw/common
	$(MAKE) -C manyclaw/ptwise/steppers/advection
	$(MAKE) -C manyclaw/ptwise/steppers/acoustics_const
	$(MAKE) -C manyclaw/ptwise/steppers/acoustics_var
	$(MAKE) -C manyclaw/ptwise/steppers/euler

lib: build
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o libmanyclaw.so $(shell find manyclaw -name "*.o")

test: lib
	$(MAKE) -C test

benchmark: lib
	$(MAKE) -C benchmark

clean:
	$(MAKE) -C manyclaw/common clean
	$(MAKE) -C manyclaw/ptwise/steppers/advection clean
	$(MAKE) -C manyclaw/ptwise/steppers/acoustics_const clean
	$(MAKE) -C manyclaw/ptwise/steppers/acoustics_var clean
	$(MAKE) -C manyclaw/ptwise/steppers/euler clean
	$(MAKE) -C test clean
	$(MAKE) -C benchmark clean
	rm -rf libmanyclaw.so
