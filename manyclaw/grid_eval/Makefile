include ../../variables.inc

OBJS=void_grid_eval.o advection_grid_eval.o acoustics_grid_eval.o euler_grid_eval.o

all: ${OBJS}

${OBJS}: ../riemann/advection_rp.h ../common/common.h ../common/data_structures.h

.PHONY:clean

clean:
	/bin/rm -rf driver *.o
