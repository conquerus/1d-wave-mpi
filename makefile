.POSIX:
CC=mpicc
CFLAGS=-std=c99 -Wall -pedantic -pthread -O3 
LDLIBS = -lm

all: main

main: main.o grid.o mpi_util.o
	$(CC) $(CFLAGS) -o ./bin/wave_1D main.o grid.o mpi_util.o $(LDFLAGS) $(LDLIBS)
main.o: main.c problem.h
grid.o: grid.c grid.h
mpi_util.o: mpi_util.c mpi_util.h

.PHONY: clean

clean:
	rm -f *.o ./bin/*

