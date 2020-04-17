#include<iostream>
#include<string>
#include<fstream>
#include<chrono>
#include<math.h>
#include<mpi.h>

extern double** block_init(double **,int, int, int*, int*,int, int);
/* Returns to each process its corresponding block
 * ==================================================================
 * Parameters : 
 * double ** : Matrix to be cut in nbProcess blocks
 * int       : Number of rows of the given matrix
 * int       : Number of columns of the given matrix
 * int*      : Pointer to the number of rows of the returned block
 * int*      : Pointer to the number of columns of the returned block
 * int       : Rank
 * int       : number of processes     
 * ==================================================================
 * Comment : The different blocks have not necessarily the same
 *           dimensions.
 */

extern double** gather_blocks(double **,int, int,int,int,int,int);
/* Rank 0 returns the gathered matrix from scattered blocks
 * The others return a NULL Pointer
 * ==================================================================
 * Parameters : 
 * double ** : Block matrix to place into the gathered matrix
 * int       : Number of rows of the given block matrix
 * int       : Number of columns of the given block matrix
 * int       : Number of rows of the gathered matrix
 * int       : Number of columns of the gathered matrix
 * int       : Rank
 * int       : number of processes     
 * ==================================================================
 * Comment : Only rank 0 will receive the gathered matrix, it is 
 *           good to know before printing it in the main program.
 */

extern double **matrix_multiplication(double **, double **, int, int, int);
/* Returns the multiplication of two matrices given as parameters.
 * ==================================================================
 * Parameters : 
 * double ** : Matrix A
 * double ** : Matrix B
 * int       : A's number of rows
 * int       : A's number of columns
 * int       : B's number of columns
 * ==================================================================
 * Comment : B's number of rows is necessarily egal to A's number of
 *           columns.
 */

extern double **parallel_matrix_multiplication(double **,double **,int,int,int,int,int);
/* Each rank returns its corresponding block of the block
 * multiplication of A and B.
 * ==================================================================
 * Parameters : 
 * double ** : Block matrix A_Block
 * double ** : Block matrix B_Block
 * int       : A_Block's number of rows
 * int       : A_Block's number of columns
 * int       : B_Block's number of columns
 * int       : Rank
 * int       : Number of Processes
 * ==================================================================
 * Comment : B_block's number of rows is necessarily egal to A's 
 *           number of columns.
 */

extern double** init_contiguous_matrix(int,int);
/* Returns a pointer to a contiguous 2D array space allocated
 * in the heap memory.
 * ==================================================================
 * int       : Number of rows
 * int       : Number of columns
 * ==================================================================
 * Comment : This function is extremely useful to pass data between
 *           processes with MPI's subroutines.
 *           Need to be cleaned with the "clean_contiguous_matrix"
 *           subroutine.
 */

extern double** init_matrix(int, int);
extern void print_matrix(double **, int, int);
extern void print_vector(double *,int);
extern void read_matrix(double **, int, int, std::string);
extern void read_vector(double *,int, std::string);
extern double* init_vector(int);
extern void clean_matrix(double **, int, int);
extern void clean_contiguous_matrix(double**,int,int);
extern void clean_vector(double *, int);
extern void matrix_copy(double **, double **, int, int);
