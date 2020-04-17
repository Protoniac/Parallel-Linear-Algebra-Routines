#include"lu.h"

void crout_lu_factorization(double ** A, double **U, double **L, int n){
    /* In coming ... */
}

void lu_factorization(double **A, double **U, double **L, int n){
    
    for(int i=0;i<n;i++){
        L[i][i] = 1;
        for(int j=0;j<n;j++)  U[i][j] = A[i][j];
    }
    for(int i=0;i<n-1;i++){
        for(int j=i+1;j<n;j++){
            L[j][i] = U[j][i]/U[i][i];
            for (int k = i;k<n;k++){
                U[j][k] -= U[i][k]*L[j][i];
            }
        }
    }
}

void parallel_lu_factorization(double **A, double **U, double **L, int n, int rank, int nbProcess){
    /* In coming ... */
}
