/*helper functions and variables for MPI*/

#ifndef MPI_UTIL
#define MPI_UTIL

#ifndef PROBLEM
#include "problem.h"
#endif

#include <mpi.h>

int ierr;
int num_procs;
int root_process;
int my_id;

MPI_Status status;

int start_mpi(int * argc, char ** argv);

/* Decide the starting and ending points in physical space for each proc. */
void allocate_to_procs(double start, double end);

/* Get the start and end values from root proc. */
void get_start_end_vals();

#endif /*MPI_UTIL*/
