CXX ?= g++
CXXFLAGS ?= 
LINK ?= $(CXX)

CXXFLAGS += -I../..
LDFLAGS = $(TBB_LDFLAGS) -L../../lib -lmanyclaw 

swirl: swirl.o
	$(LINK) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

clean:
	-rm -rf swirl swirl.o

clobber: clean
	$(MAKE) -C ../../ clean

new: clobber
	$(MAKE) -C ../../
	$(MAKE) swirl
