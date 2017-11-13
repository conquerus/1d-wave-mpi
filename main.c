#ifndef PROBLEM
#include "problem.h"
#endif

#ifndef MPI_UTIL
#include "mpi_util.h"
#endif

#ifndef GRID
#include "grid.h"
#endif

int main(int argc, char* argv[])
{
  /* MPI */
  extern int my_id;
  extern int ierr;
  extern int num_procs;
  extern int root_process;

  /* Time */
  double time = INIT_TIME;

  /* Space */
  double start = 0;
  double end = LENGTH;
  extern int N_local;  

  start_mpi(&argc, argv);

  N_local = N/num_procs;
  
  if (my_id == root_process) {
    allocate_to_procs(start, end);
  }

  /* each proc gests its physical space allocations */
  get_start_end_vals(&start, &end);
  
  /* Grid */
  struct point ghost_left;
  struct point point_array[N_local];
  struct point ghost_right;
  
  /* initialize */
  init_array(point_array, start, end);

  /* main time marching loop */
  while (time <= MAX_TIME) {
    apply_BCs(point_array, &ghost_left, &ghost_right, my_id);
    update_grid(point_array, &ghost_left, &ghost_right);
    update_storage(point_array);
    
    time += DT;
  }

  output(point_array);
  
  ierr = MPI_Finalize();
  return 0;
}
