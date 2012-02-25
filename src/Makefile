include variables.inc
include rules.inc

all:
	$(MAKE) -C common
	$(MAKE) -C ptwise/steppers/advection
	$(MAKE) -C drivers

clean:
	$(MAKE) -C common clean
	$(MAKE) -C ptwise/steppers/advection clean
	$(MAKE) -C drivers clean
