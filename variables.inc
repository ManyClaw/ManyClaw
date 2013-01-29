CXX=g++
OMPCFLAGS=-fopenmp
CXXFLAGS=-I ./ -O3 -Wall $(OMPCFLAGS)
LDFLAGS= -ltbb


