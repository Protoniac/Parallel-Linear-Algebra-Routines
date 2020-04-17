#include"lu.h"
#include"matrix.h"

int main(int argc, char *argv[]){
    
    int nrowA = 1000;
    int ncolA = 1500;
    int ncolB = 765;
    
    double **A;
    double **B;
    
    double **AB;
    
    int rank, nbProcess;
    
    MPI_Init(&argc,&argv);
    
    MPI_Comm_size(MPI_COMM_WORLD,&nbProcess);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    if(rank == 0){
        
        A = init_contiguous_matrix(nrowA,ncolA);
        B = init_contiguous_matrix(ncolA,ncolB);
        
        read_matrix(A,nrowA,ncolA,"matrix_samples/1000x1500_1.txt");
        read_matrix(B,ncolA,ncolB,"matrix_samples/1500x765_1.txt");
    }
    else{
        A = (double **)malloc(sizeof(double**));
        B = (double **)malloc(sizeof(double**));
    }
    
    auto start = std::chrono::steady_clock::now();
    
    int nrowBlockA, ncolBlockA, ncolBlockB;
    
    double **block_A = block_init(A,nrowA,ncolA,&nrowBlockA,&ncolBlockA,rank,nbProcess);
    
    double **block_B = block_init(B,ncolA,ncolB,&ncolBlockA,&ncolBlockB,rank,nbProcess);
    
    double **AB_block = parallel_matrix_multiplication(block_A,block_B,nrowBlockA,ncolBlockA,ncolBlockB,rank,nbProcess);
    
    double **AB_gathered = gather_blocks(AB_block,nrowBlockA,ncolBlockB,nrowA,ncolB,rank,nbProcess);
    
    auto end = std::chrono::steady_clock::now();
    
    if(rank == 0){
        
//         print_matrix(AB_gathered,nrowA,ncolB);
        AB = matrix_multiplication(A,B,nrowA,ncolA,ncolB);
//         print_matrix(AB,nrowA,ncolB);
        
        std::cout << "Elapsed Time : "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()
        << " ms. " << std::endl;
        
        clean_contiguous_matrix(A,nrowA,ncolA);
        clean_contiguous_matrix(B,ncolA,ncolB);
        clean_contiguous_matrix(AB,nrowA,ncolB);
        clean_contiguous_matrix(AB_gathered,nrowA,ncolB);
        
    }
    else{
        free(A);
        free(B);
    }
    
    clean_contiguous_matrix(block_A,nrowBlockA,ncolBlockA);
    clean_contiguous_matrix(block_B,ncolBlockA,ncolBlockB);
    clean_contiguous_matrix(AB_block,nrowBlockA,ncolBlockB);
    
    MPI_Finalize();
 
    return 0;
}
