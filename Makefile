include variables.inc

all: lib test benchmark

lib:
	$(MAKE) -C src/common
	$(MAKE) -C src/ptwise/steppers/advection
	$(MAKE) -C src/ptwise/steppers/acoustics_const
	$(MAKE) -C src/ptwise/steppers/acoustics_var
	$(MAKE) -C src/ptwise/steppers/euler

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
