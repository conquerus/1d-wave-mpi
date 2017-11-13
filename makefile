.POSIX:
CC = mpicc
CFLAGS = -std=c99 -Wall -pedantic -pthread -O3 
LDFLAGS = -lm
LDLIBS = 

OBJECTS = main.o grid.o mpi_util.o

all: main

main: $(OBJECTS)
	$(CC) $(CFLAGS) -o ./bin/wave_1D $(OBJECTS) $(LDFLAGS) $(LDLIBS)

main.o: main.c problem.h
grid.o: grid.c grid.h
mpi_util.o: mpi_util.c mpi_util.h

.PHONY: clean

clean:
	rm -f *.o ./bin/*

