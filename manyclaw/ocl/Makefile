MANYCLAW_DIR = ../../
include $(MANYCLAW_DIR)/variables.inc

OCL_SRC = CLManager.cpp
OCL_OBJ = $(subst .cpp,.o, $(OCL_SRC))


all: $(OCL_OBJ)

CLManager.o: CLManager.h

clean:
	rm -rf common.o data_structures.o
