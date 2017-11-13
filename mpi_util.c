#include "mpi_util.h"

int start_mpi(int * argc, char ** argv) {
  ierr = MPI_Init(argc, &argv);

  /* find out my process ID, and how many processes were started. */
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  root_process = 0;
  
  return ierr;
}

/* Decide the starting and ending points in physical space for each proc. */
void allocate_to_procs(double start, double end) {
  for (int id = 0; id < num_procs; id++) {
    start = (id)*(LENGTH/num_procs);
    end   = (id+1)*(LENGTH/num_procs);

    /* send this info to the children */
    ierr = MPI_Send(&start, 1, MPI_DOUBLE, id, 1, MPI_COMM_WORLD);
    ierr = MPI_Send(&end,   1, MPI_DOUBLE, id, 1, MPI_COMM_WORLD);
  }
}

/* Get the start and end values from root proc. */
void get_start_end_vals(double * start, double * end) {
  ierr = MPI_Recv(start, 1, MPI_DOUBLE, root_process, 1, MPI_COMM_WORLD, &status);
  ierr = MPI_Recv(end,   1, MPI_DOUBLE, root_process, 1, MPI_COMM_WORLD, &status);
}

