include ../../../variables.inc

OBJS = euler_rp_step_serial.o  \
       euler_rp_step_serial_cellwise_functor.o \
       euler_rp_step_serial_cellwise.o \
       euler_rp_step_serial_tiled.o \
       euler_rp_step_omp.o \
       euler_rp_step_omp2.o \
       euler_rp_step_tbb.o

all: $(OBJS)

$(OBJS): ../../riemann/euler_rp.h ../../common/common.h ../../common/data_structures.h

.PHONY: clean

clean:
	/bin/rm -rf $(OBJS)
