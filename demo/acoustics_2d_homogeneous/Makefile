CXX ?= g++
CXXFLAGS ?= 
LINK ?= $(CXX)

CXXFLAGS += -I../..
LDFLAGS = $(TBB_LDFLAGS) -L../../lib -lmanyclaw 

acoustics_homogeneous: acoustics_homogeneous.o
	$(LINK) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

clean:
	-rm -rf acoustics_homogeneous acoustics_homogeneous.cpp

clobber: clean
	$(MAKE) -C ../../ clean
