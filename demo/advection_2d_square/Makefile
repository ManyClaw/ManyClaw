CXX ?= g++
CXXFLAGS ?= 
LINK ?= $(CXX)

CXXFLAGS += -I../..
LDFLAGS = $(TBB_LDFLAGS) -L../../lib -lmanyclaw 

advection: advection.o
	$(LINK) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

clean:
	-rm -rf advection advection.o

clobber: clean
	$(MAKE) -C ../../ clean
