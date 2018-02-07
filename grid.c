#include "grid.h"

unsigned int N_local = 0;

/* Make point array and assign initial condition */
void init_array(struct point point_array[], double start, double stop)
{
  double step = (stop-start)/N_local;
  double x = start;
  
  for (unsigned int i = 0; i <= N_local; i++) {
    point_array[i].x = x;
    point_array[i].val = INITIAL_CONDITION(x);
    point_array[i].val_old = point_array[i].val;
    x += step;
  }
  return;
}

/* Apply the BCs and communicate between procs. */
void apply_BCs(struct point point_array[],
               struct point* ghost_left,
               struct point* ghost_right,
               int my_id)
{
  if (my_id == 0) {
    /* send ghost cells */
    ierr = MPI_Send(&(point_array[N_local-1].val), 1, MPI_DOUBLE, my_id+1, 1, MPI_COMM_WORLD);
    /* receive ghost cells */
    ierr = MPI_Recv(&(ghost_right->val), 1, MPI_DOUBLE, my_id+1, 1, MPI_COMM_WORLD, &status);
    /* BC */
    (*ghost_left).x = -DX;
    (*ghost_left).val = 0;
    (*ghost_left).val_old = 0;
    (*ghost_left).val_new = 0;
  }
  else if (my_id == num_procs-1) {
    /* send ghost cells */
    ierr = MPI_Send(&(point_array[0].val),         1, MPI_DOUBLE, my_id-1, 1, MPI_COMM_WORLD);
    /* receive ghost cells */
    ierr = MPI_Recv(&(ghost_left->val),  1, MPI_DOUBLE, my_id-1, 1, MPI_COMM_WORLD, &status);
    /* BC */
    (*ghost_right).x = LENGTH+DX;
    (*ghost_right).val = 0;
    (*ghost_right).val_old = 0;
    (*ghost_right).val_new = 0;
  }
  else {
    /* send ghost cells */
    ierr = MPI_Send(&(point_array[0].val),         1, MPI_DOUBLE, my_id-1, 1, MPI_COMM_WORLD);
    ierr = MPI_Send(&(point_array[N_local-1].val), 1, MPI_DOUBLE, my_id+1, 1, MPI_COMM_WORLD);
    /* receive ghost cells */
    ierr = MPI_Recv(&(ghost_right->val), 1, MPI_DOUBLE, my_id+1, 1, MPI_COMM_WORLD, &status);
    ierr = MPI_Recv(&(ghost_left->val),  1, MPI_DOUBLE, my_id-1, 1, MPI_COMM_WORLD, &status);
  }
}

/* Perform a single time step using first-order spacial discretization. */
void update_grid(struct point point_array[], struct point* ghost_left, struct point* ghost_right) {
  double CFL2 = ((SPEED*DT)/DX)*((SPEED*DT)/DX);
  
  point_array[0].val_new = CFL2*(point_array[1].val + ghost_left->val) + 2.0*(1.0-CFL2)*point_array[0].val - point_array[0].val_old;
  
  for (unsigned int i = 1; i <= (N_local - 2); i++) {
    point_array[i].val_new = CFL2*(point_array[i+1].val + point_array[i-1].val) + 2.0*(1.0-CFL2)*point_array[i].val - point_array[i].val_old;
  }
  
  point_array[N_local-1].val_new = CFL2 * (ghost_right->val + point_array[N_local-2].val) + 2.0*(1.0-CFL2)*point_array[N_local-1].val - point_array[N_local-1].val_old;
}

/*update storage  val_old = val -> val = val_new  */
void update_storage(struct point point_array[]) {
  for (unsigned int i = 0; i < N_local; i++) {
    point_array[i].val_old = point_array[i].val;
    point_array[i].val = point_array[i].val_new;
  }
}

void output(struct point point_array[]) {
  for (unsigned int i = 0; i < N_local; i++) {
    printf("%f; %f\n", point_array[i].x, point_array[i].val);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
