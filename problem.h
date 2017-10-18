/* The problem definition */

#ifndef PROBLEM
#define PROBLEM

#include <math.h>

#define PI 3.1415926

#define LENGTH 6.0
#define DX 0.001001001
#define N (unsigned int)(LENGTH/DX + 1.0)

#define INIT_TIME 0.0
#define MAX_TIME 30.0
#define DT 0.0001

#define SPEED 4.0

#define INITIAL_CONDITION(x) (double)(0.5 - 0.5*cos((2.0*PI/LENGTH)*(x)))

#endif /*PROBLEM*/