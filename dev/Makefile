# Makefile for Dispatcher/dev

#Makefile for Dispatcher/dev folder

TARGET = bdev
PROCS = 2
CFLAGS = -O0
DEBUG = -g

$(TARGET) : $(TARGET).cxx
	mpic++ -std=c++17 $(TARGET).cxx -o ../bin/$(TARGET)
	
millerrabin.o : millerrabin.cxx
	mpic++ -c millerrabin.cxx -o millerrabin.o

PHONY : run

run :
	mpirun -np $(PROCS) -H localhost:12 ../bin/$(TARGET)
