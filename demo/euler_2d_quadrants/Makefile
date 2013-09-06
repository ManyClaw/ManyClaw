CXX ?= g++
CXXFLAGS ?= 
LINK ?= $(CXX)

CXXFLAGS += -I../.. -g
LDFLAGS = $(TBB_LDFLAGS) -L../../lib -lmanyclaw 

quadrants: quadrants.o
	$(LINK) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

clean:
	-rm -rf quadrants quadrants.o

clobber: clean
	$(MAKE) -C ../../ clean

new: clobber
	$(MAKE) -C ../../
	$(MAKE) quadrants
