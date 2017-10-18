#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>

#define PI 3.1415926

#define LENGTH 6.0
#define DX 0.001001001
//#define DX 0.005005005

#define MAX_TIME 30.0
#define DT 0.0001

#define SPEED 4.0

struct point {
  double x;
  double val_old;
  double val;
  double val_new;
};

double initial_condition(double x)
{
  return 0.5 - 0.5*cos((2.0*PI/LENGTH)*x);
}

void init_array_step(struct point p_point[], int start, double step, unsigned int N)
{
  double x = start;
  unsigned int i;

  for (i = 0; i <= N; i++) {
    p_point[i].x = x;
    p_point[i].val = initial_condition(x);
    p_point[i].val_old = p_point[i].val;
    x += step;
  }
  return;
}

void init_array_end_points(struct point p_point[], double start, double stop, unsigned int N)
{
  double x = start;
  double step = (stop-start)/N;
    
  unsigned int i;
  for (i = 0; i <= N; i++) {
    p_point[i].x = x;
    p_point[i].val = initial_condition(x);
    p_point[i].val_old = p_point[i].val;
    x += step;
  }
  return;
}

int main(int argc, char** argv)
{
  double time = 0.0;
  unsigned int N = LENGTH/DX + 1;

  double CFL2 = ((SPEED*DT)/DX)*((SPEED*DT)/DX);
    

  int i;

  volatile int inf=0;
  
  double start = 0;
  double end = LENGTH;
  
  int ierr, num_procs, my_id;
  int root_process = 0;

  MPI_Status status;
  ierr = MPI_Init(&argc, &argv);

  /* find out MY process ID, and how many processes were started. */

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  //while (inf==0); //wait for the debugger
  //MPI_Barrier(MPI_COMM_WORLD);
  
  if (my_id == root_process) {
    printf("root process hi");

    for (int id = 0; id < num_procs; id++) {
      start = (id)*(LENGTH/num_procs);
      end   = (id+1)*(LENGTH/num_procs);

      //send this info to the children
      ierr = MPI_Send(&start, 1, MPI_DOUBLE, id, 1, MPI_COMM_WORLD);
      ierr = MPI_Send(&end,   1, MPI_DOUBLE, id, 1, MPI_COMM_WORLD);
    }
  }

  //MPI_Barrier(MPI_COMM_WORLD);
  ierr = MPI_Recv(&start, 1, MPI_DOUBLE, root_process, 1, MPI_COMM_WORLD, &status);
  ierr = MPI_Recv(&end,   1, MPI_DOUBLE, root_process, 1, MPI_COMM_WORLD, &status);

  unsigned int N_local = N/num_procs;
  struct point ghost_left;
  struct point point_array[N_local];
  struct point ghost_right;
  
  printf("process %i: start %f end %f size %d\n", my_id, start, end, N_local);
  
  /*initialize*/
  init_array_end_points(point_array, start, end, N_local);

  while (time <= MAX_TIME) {

    if (my_id == 0) {
      //send ghost cells
      ierr = MPI_Send(&(point_array[N_local-1].val), 1, MPI_DOUBLE, my_id+1, 1, MPI_COMM_WORLD);

      //receive ghost cells
      ierr = MPI_Recv(&(ghost_right.val), 1, MPI_DOUBLE, my_id+1, 1, MPI_COMM_WORLD, &status);
      ghost_left.x = -DX;
      ghost_left.val = 0;
      ghost_left.val_old = 0;
      ghost_left.val_new = 0;
    }
    else if (my_id == num_procs-1) {
      //send ghost cells
      ierr = MPI_Send(&(point_array[0].val),         1, MPI_DOUBLE, my_id-1, 1, MPI_COMM_WORLD);

      //receive ghost cells
      ghost_right.x = LENGTH+DX;
      ghost_right.val = 0;
      ghost_right.val_old = 0;
      ghost_right.val_new = 0;
      ierr = MPI_Recv(&(ghost_left.val),  1, MPI_DOUBLE, my_id-1, 1, MPI_COMM_WORLD, &status);
    }
    else {
      //send ghost cells
      ierr = MPI_Send(&(point_array[0].val),         1, MPI_DOUBLE, my_id-1, 1, MPI_COMM_WORLD);
      ierr = MPI_Send(&(point_array[N_local-1].val), 1, MPI_DOUBLE, my_id+1, 1, MPI_COMM_WORLD);

      //receive ghost cells
      ierr = MPI_Recv(&(ghost_right.val), 1, MPI_DOUBLE, my_id+1, 1, MPI_COMM_WORLD, &status);
      ierr = MPI_Recv(&(ghost_left.val),  1, MPI_DOUBLE, my_id-1, 1, MPI_COMM_WORLD, &status);
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    
    point_array[0].val_new = CFL2*(point_array[1].val + ghost_left.val) + 2.0*(1.0-CFL2)*point_array[0].val - point_array[0].val_old;
        
    for (i = 1; i <= (N_local - 2); i++) {
      point_array[i].val_new = CFL2*(point_array[i+1].val + point_array[i-1].val) + 2.0*(1.0-CFL2)*point_array[i].val - point_array[i].val_old;
    }
    
    point_array[N_local-1].val_new = CFL2*(ghost_right.val + point_array[N_local-2].val) + 2.0*(1.0-CFL2)*point_array[N_local-1].val - point_array[N_local-1].val_old;
    
    for (i = 0; i <= (N_local - 1); i++) {
      /*update storage  val_old = val  val = val_new  */
      point_array[i].val_old = point_array[i].val;
      point_array[i].val = point_array[i].val_new;
    }

    time += DT;
    //MPI_Barrier(MPI_COMM_WORLD);
  }

  //output
  for (i = 0; i < (N_local - 1); i++) {
    printf("%f, %f\n", point_array[i].x, point_array[i].val);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  ierr = MPI_Finalize();
  return 0;
}
