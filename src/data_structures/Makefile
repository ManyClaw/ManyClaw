include ../variables.inc

all: grid.o

grid.o: grid.cpp grid.h
	$(CXX) $(CXXFLAGS) -o $@ -c grid.cpp

clean:
	rm -rf grid.o
