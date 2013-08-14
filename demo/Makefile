EXAMPLES = \
	advection \
	swirl \
	euler_quadrants

all: $(EXAMPLES)

advection:
	$(MAKE) -C advection_2d_square

swirl:
	$(MAKE) -C advection_2d_swirl

euler:
	$(MAKE) -C euler_2d_quadrants

clean:
	$(MAKE) -C advection_2d_square clean
	$(MAKE) -C advection_2d_swirl clean
	$(MAKE) -C euler_2d_quadrants clean

clobber: clean
	$(MAKE) -C ../ clean
