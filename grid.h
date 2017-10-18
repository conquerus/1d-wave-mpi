#ifndef GRID
#define GRID

#include <stdio.h>

#ifndef PROBLEM
#include "problem.h"
#endif

#ifndef MPI_UTIL
#include "mpi_util.h"
#endif

int N_local;

struct point {
  double x;
  double val_old;
  double val;
  double val_new;
};

void init_array(struct point p_point[], double start, double stop);

/* Apply the BCs and communicate between procs. */
void apply_BCs(struct point point_array[],
               struct point* ghost_left,
               struct point* ghost_right,
               int my_id);

/* Perform a single time step using first order spacial discretiztion. */
void update_grid(struct point point_array[],
                 struct point* ghost_left,
                 struct point* ghost_right);

/*update storage  val_old = val  val = val_new  */
void update_storage(struct point point_array[]);

void output(struct point point_array[]);

#endif /*GRID*/
