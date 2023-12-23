# Makefile for Dispatcher/dev

TARGET = adev
PROCS = 12

$(TARGET) : adev.cxx
	mpic++ -std=c++17 -g adev.cxx -o ../bin/$(TARGET)

PHONY : run

run :
	mpirun -np $(PROCS) -H localhost:12 ../bin/$(TARGET)
