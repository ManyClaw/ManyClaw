include ../../variables.inc

OBJS=advection_grid_eval.o

all: ${OBJS}

${OBJS}: ../riemann/advection_rp.h ../common/common.h ../common/data_structures.h

.PHONY:clean

clean:
	/bin/rm -rf driver *.o
