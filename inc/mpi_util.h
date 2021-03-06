/*helper functions and variables for MPI*/

#ifndef MPI_UTIL
#define MPI_UTIL

#ifndef PROBLEM
#include "problem.h"
#endif

#include <mpi.h>

extern int ierr;
extern int num_procs;
extern int root_process;
extern int my_id;

MPI_Status status;

int start_mpi(int * argc, char ** argv);

/* Decide the starting and ending points in physical space for each proc. */
void allocate_to_procs(double start, double end);

/* Get the start and end values from root proc. */
void get_start_end_vals(double * start, double * end);

#endif /*MPI_UTIL*/
