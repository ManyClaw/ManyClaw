include variables.inc
include rules.inc

all:
	$(MAKE) -C common
	$(MAKE) -C ptwise/steppers/advection
	$(MAKE) -C ptwise/steppers/acoustics_const
	$(MAKE) -C ptwise/steppers/euler
	$(MAKE) -C drivers

clean:
	$(MAKE) -C common clean
	$(MAKE) -C ptwise/steppers/advection clean
	$(MAKE) -C ptwise/steppers/acoustics_const clean
	$(MAKE) -C ptwise/steppers/euler clean
	$(MAKE) -C drivers clean
