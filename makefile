.POSIX:
IDIR = ./inc
CC = mpicc
CFLAGS = -std=c99 -Wall -Wextra -pedantic -pthread -O3 -march=native
CFLAGS += -I$(IDIR)
LDFLAGS = -lm
LDLIBS = 

OBJECTS = main.o grid.o mpi_util.o

all: main

main: $(OBJECTS)
	$(CC) $(CFLAGS) -o ./bin/wave_1D $(OBJECTS) $(LDFLAGS) $(LDLIBS)

main.o: main.c ./inc/problem.h
grid.o: grid.c ./inc/grid.h
mpi_util.o: mpi_util.c ./inc/mpi_util.h

.PHONY: clean

clean:
	rm -f *.o ./bin/wave_1D

