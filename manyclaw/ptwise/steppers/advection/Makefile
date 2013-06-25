include ../../../../variables.inc

OBJS=advection_rp_step_serial.o advection_rp_step_serial_cellwise.o advection_rp_step_serial_tiled.o advection_rp_step_omp.o advection_rp_step_tbb.o advection_rp_step_serial_cellwise_indexer.o

all: ${OBJS}

${OBJS}: ../../riemann/advection_rp.h ../../../common/common.h ../../../common/data_structures.h

.PHONY:clean

clean:
	/bin/rm -rf driver *.o
