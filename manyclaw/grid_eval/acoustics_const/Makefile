include ../../../variables.inc

OBJS = acoustics_const_rp_step_serial.o \
       acoustics_const_rp_step_serial_cellwise.o \
       acoustics_const_rp_step_serial_tiled.o \
       acoustics_const_rp_step_omp.o \
       acoustics_const_rp_step_tbb.o

all: $(OBJS)

$(OBJS): ../../riemann/acoustics_const_rp.h ../../common/common.h ../../common/data_structures.h

.PHONY: clean

clean:
	/bin/rm -rf $(OBJS)
