#include"matrix.h"

double **block_init(double **M, int nrow, int ncol, int *pnrowBlock, int *pncolBlock, int rank, int nbProcess){
    
    int subsizes[2];
    int starts[2] = {0,0};
    
    int sendcounts[nbProcess];
    int senddispls[nbProcess];
    
    int recvcounts[nbProcess];
    int recvdispls[nbProcess];
    
    MPI_Datatype recvtypes[nbProcess];
    
    int globalsizes[2] = {nrow,ncol};
    
    int total_shift = 0;
    int row_shift = 0;
    
    int row_remaining = nrow%(int)sqrt(nbProcess);
    int col_remaining = ncol%(int)sqrt(nbProcess);
    
    int row_part = (nrow-row_remaining)/(int)sqrt(nbProcess);
    int col_part = (ncol-col_remaining)/(int)sqrt(nbProcess);
    
    double **block;
    
    MPI_Datatype blocktypes[nbProcess];
    
    for(int i=0;i<(int)sqrt(nbProcess);i++){
        for(int j=0;j<(int)sqrt(nbProcess);j++){
            subsizes[0] = row_part;
            subsizes[1] = col_part;
            if(row_remaining > 0) subsizes[0] +=1;
            if(col_remaining > 0){
                subsizes[1] +=1;
                col_remaining -= 1;
            }
            
            MPI_Type_create_subarray(2,globalsizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&blocktypes[i*(int)sqrt(nbProcess)+j]);
            MPI_Type_commit(&blocktypes[i*(int)sqrt(nbProcess)+j]);
            
            if(rank == i*(int)sqrt(nbProcess)+j){
                
                *pnrowBlock = subsizes[0];
                *pncolBlock = subsizes[1];
                block = init_contiguous_matrix(*pnrowBlock,*pncolBlock);
                
            }
            
            if(rank == 0){
//                 std::cout << subsizes[0] << "," << subsizes[1] << std::endl;
                sendcounts[i*(int)sqrt(nbProcess)+j] = 1;
                senddispls[i*(int)sqrt(nbProcess)+j] = (total_shift)*sizeof(double);
                total_shift += subsizes[1];
                row_shift += subsizes[1];
            }
            else{
                sendcounts[i*(int)sqrt(nbProcess)+j] = 0;
                senddispls[i*(int)sqrt(nbProcess)+j] = 0;
            }

            recvdispls[i*(int)sqrt(nbProcess)+j] = 0;
            recvcounts[i*(int)sqrt(nbProcess)+j] = 0;
            recvtypes[i*(int)sqrt(nbProcess)+j] = MPI_DOUBLE;
            
        }
        if(rank==0){
            
            total_shift += row_shift*(subsizes[0]-1);
            row_shift = 0;
        }
        row_remaining -= 1;
        col_remaining = ncol%(int)sqrt(nbProcess);
    }
    recvcounts[0] = (*pnrowBlock)*(*pncolBlock);
    
    MPI_Alltoallw(&(M[0][0]),sendcounts,senddispls,blocktypes,&(block[0][0]),recvcounts,recvdispls,recvtypes,MPI_COMM_WORLD);

    return block;
}

double ** matrix_multiplication(double **A, double **B, int nrowA, int nrowB, int ncolB){
    double **AB = init_contiguous_matrix(nrowA,ncolB);
    for(int i = 0;i<nrowA;i++){
        for(int j = 0;j<ncolB;j++){
            AB[i][j] = 0;
            for(int k = 0;k<nrowB;k++){
                AB[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return AB;
}

void matrix_copy(double **B, double **A, int nrow, int ncol){
    for(int i = 0;i<nrow;i++){
        for(int j = 0;j<ncol;j++){
            B[i][j] = A[i][j];
        }
    }
}

double ** gather_blocks(double **AB_block, int nrowBlockA, int ncolBlockB,int nrowA, int ncolB, int rank, int nbProcess){
    
    double **AB_Gathered;

    int globalsizes[2] = {nrowA,ncolB};
    int starts[2] = {0,0};
    int subsizes[2];
    
    int total_shift = 0;
    int row_shift = 0;
    
    int sendcounts[nbProcess];
    int senddispls[nbProcess];
    MPI_Datatype sendtypes[nbProcess];
    
    MPI_Datatype blocktypes[nbProcess];
    int recvcounts[nbProcess];
    int recvdispls[nbProcess];
    
    for(int i = 0;i<(int)sqrt(nbProcess);i++){
        for(int j = 0;j<(int)sqrt(nbProcess);j++){
            
            if(rank == i*(int)sqrt(nbProcess) + j){
                subsizes[0] = nrowBlockA;
                subsizes[1] = ncolBlockB;
            }
            
            MPI_Bcast(subsizes,2,MPI_INT,i*(int)sqrt(nbProcess) + j,MPI_COMM_WORLD);
            
            MPI_Type_create_subarray(2,globalsizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&blocktypes[i*(int)sqrt(nbProcess)+j]);
            
            MPI_Type_commit(&blocktypes[i*(int)sqrt(nbProcess)+j]);
            
            if(rank == 0) {
                recvcounts[i*(int)sqrt(nbProcess)+j] = 1;
                recvdispls[i*(int)sqrt(nbProcess)+j] = total_shift*sizeof(double);
                total_shift += subsizes[1];
                row_shift += subsizes[1];
            }
            else{
                recvcounts[i*(int)sqrt(nbProcess)+j] = 0;
                recvdispls[i*(int)sqrt(nbProcess)+j] = 0;
            }
            sendcounts[i*(int)sqrt(nbProcess) +j] = 0;
            senddispls[i*(int)sqrt(nbProcess) +j] = 0;
            sendtypes[i*(int)sqrt(nbProcess) +j] = MPI_DOUBLE;
        }
        total_shift += row_shift*(subsizes[0]-1);
        row_shift = 0;
    }
    
    sendcounts[0] = nrowBlockA*ncolBlockB;
    
    if(rank ==0){
        AB_Gathered = init_contiguous_matrix(nrowA,ncolB);
    }
    MPI_Alltoallw(&(AB_block[0][0]),sendcounts,senddispls,sendtypes,&(AB_Gathered[0][0]),recvcounts,recvdispls,blocktypes,MPI_COMM_WORLD);
    
    if(rank == 0) return AB_Gathered;
    else return (double **)0;
    
}

double ** parallel_matrix_multiplication(double **A_block,double **B_block,int nrowBlockA,int ncolBlockA,int ncolBlockB,int rank,int nbProcess){

    int I = (int)(rank/(int)sqrt(nbProcess)) + 1;
    int J = (rank%(int)sqrt(nbProcess)) + 1;
    
    MPI_Comm COMM_ROW;
    MPI_Comm COMM_COL;
    
    MPI_Comm_split(MPI_COMM_WORLD,J,I,&COMM_COL);
    MPI_Comm_split(MPI_COMM_WORLD,I,J,&COMM_ROW);
    
    int ncolBlockAnrowBlockB_temp;
    
    double **A_block_temp;
    double **B_block_temp;
    
    double **AB_block = init_contiguous_matrix(nrowBlockA,ncolBlockB);
        
    for(int k = 1;k<(int)sqrt(nbProcess)+1;k++){
        
        if(k==I){
            
            B_block_temp = init_contiguous_matrix(ncolBlockA,ncolBlockB);
            matrix_copy(B_block_temp,B_block,ncolBlockA,ncolBlockB);
            ncolBlockAnrowBlockB_temp = ncolBlockA;
            
        }
        
        MPI_Bcast(&ncolBlockAnrowBlockB_temp,1,MPI_INT,k-1,COMM_COL);
        
        if(k!=I) B_block_temp = init_contiguous_matrix(ncolBlockAnrowBlockB_temp,ncolBlockB);
        
        MPI_Bcast(&(B_block_temp[0][0]),ncolBlockAnrowBlockB_temp*ncolBlockB,MPI_DOUBLE,k-1,COMM_COL);
        
        if(k==J){
            
            A_block_temp = init_contiguous_matrix(nrowBlockA,ncolBlockA);
            matrix_copy(A_block_temp,A_block,nrowBlockA,ncolBlockA); //Idem
            ncolBlockAnrowBlockB_temp = ncolBlockA;
            
        }
        
        MPI_Bcast(&ncolBlockAnrowBlockB_temp,1,MPI_INT,k-1,COMM_ROW);
        
        if(k!=J) A_block_temp = init_contiguous_matrix(nrowBlockA,ncolBlockAnrowBlockB_temp);
        
        MPI_Bcast(&(A_block_temp[0][0]),nrowBlockA*ncolBlockAnrowBlockB_temp,MPI_DOUBLE,k-1,COMM_ROW);
        
        for(int i = 0;i<nrowBlockA;i++){
                for(int j = 0;j<ncolBlockB;j++){
                    for(int l = 0;l<ncolBlockAnrowBlockB_temp;l++){
                        AB_block[i][j] += A_block_temp[i][l]*B_block_temp[l][j];
                    }
                }
            }
        clean_contiguous_matrix(A_block_temp,nrowBlockA,ncolBlockAnrowBlockB_temp);
        clean_contiguous_matrix(B_block_temp,ncolBlockAnrowBlockB_temp,ncolBlockB);
    }
    
    return AB_block;
}

double* init_vector(int n){
    /*Space allocated here will be automatically contiguous */
    return (double *)malloc(sizeof(double)*n);
}

double** init_matrix(int nrow, int ncol){
    double ** matrix = (double **)malloc(sizeof(double *)*nrow);
    for(int i = 0;i<nrow;i++) matrix[i] = init_vector(ncol);
    return matrix;
}

double** init_contiguous_matrix(int nrow,int ncol){
    /* Allocation of the total space required in a contiguous vector */
    double *space = (double*)malloc(sizeof(double) * nrow * ncol);
    /* Allocating the space for nrow pointers*/
    double **matrix = (double **)malloc(sizeof(double*) * nrow);
    /* Linking each row of the matrix to the contiguous space */
    for(int i=0;i<nrow;i++) matrix[i] = &(space[i*ncol]);
    return matrix;
    
}

void read_vector(double *v, int n, std::string file){
    std::ifstream data;
    std::string line;
    data.open(file);
    if(data.is_open()){
        int i = 0;
        while(std::getline(data,line) && i<n){
            v[i] = stod(line);
            i++;
        }
    }
    data.close();
}

void read_matrix(double **m, int nrow, int ncol, std::string file){
    std::ifstream data;
    std::string line;
    std::string delimiter = ",";
    size_t pos = 0;
    int j;
    data.open(file);
    if(data.is_open()){
        int i = 0;
        while(std::getline(data,line) && i<nrow){
            j=0;
            while((pos = line.find(delimiter)) != std::string::npos) {
                m[i][j] = stod(line.substr(0,pos));
                line.erase(0,pos+delimiter.length());
                j++;
            }
            m[i][j] = stod(line);
            i++;
        }
    }
    data.close();
}

void print_vector(double *v, int n){
    std::cout << "Vector of size " << n << " : \n[";
    for(int i = 0; i<n; i++) std::cout << " " << v[i] << " ";
    std::cout<< "]" << std::endl;
}

void print_matrix(double **m, int nrow, int ncol){
    std::cout << "Matrix " << nrow << "x" << ncol << " : " << std::endl;
    for(int i = 0; i<nrow; i++){
        std::cout << "[";
        for(int j = 0; j<ncol; j++){
            std::cout << " "<< m[i][j] << " ";
        }
        std::cout << "]" << std::endl;
    }
}

void clean_vector(double *v, int n)
{
    free(v);
}

void clean_matrix(double **m, int nrow, int ncol){
    for(int i = 0; i < nrow;i++) clean_vector(m[i],ncol);
    free(m);
}

void clean_contiguous_matrix(double **m, int nrow, int ncol){
    free(&m[0][0]);
    free(m);
}
