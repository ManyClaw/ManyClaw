include ../../variables.inc

COMMON_SRC = common.cpp data_structures.cpp
COMMON_OBJ = $(subst .cpp,.o, $(COMMON_SRC))

all: $(COMMON_OBJ)

common.o: common.h

data_structures.o: data_structures.h common.h

clean:
	rm -rf common.o data_structures.o
