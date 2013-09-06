CXX ?= g++
CXXFLAGS ?= 
LINK ?= $(CXX)

CXXFLAGS += -I../..
LDFLAGS = $(TBB_LDFLAGS) -L../../lib -lmanyclaw 

acoustics: acoustics.o
	$(LINK) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

clean:
	-rm -rf acoustics acoustics.cpp

clobber: clean
	$(MAKE) -C ../../ clean

new: clobber
	$(MAKE) -C ../../
	$(MAKE) acoustics

