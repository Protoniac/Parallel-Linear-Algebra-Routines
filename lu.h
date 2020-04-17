#include"matrix.h"
#include<mpi.h>

extern void crout_lu_factorization(double **,double **, double **,int);
extern void lu_factorization(double **,double **,double **,int);
extern void parallel_lu_factorization(double **,double **, double **, int, int, int);
