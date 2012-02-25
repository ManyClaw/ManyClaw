CXX=g++
CXXFLAGS=-I./ -O3 -Wall
LDFLAGS= -ltbb
ISPC=ispc
ISPCFLAGS=-O2 --arch=x86-64 --target=sse2 #,sse4-x2,avx-x2
OMPCFLAGS=-fopenmp
