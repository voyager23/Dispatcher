# Makefile for Dispatcher/dev

#Makefile for Dispatcher/dev folder

TARGET = bdev
PROCS = 12

$(TARGET) : $(TARGET).cxx
	mpic++ -std=c++17 -g $(TARGET).cxx -o ../bin/$(TARGET)

PHONY : run

run :
	mpirun -np $(PROCS) -H localhost:12 ../bin/$(TARGET)
