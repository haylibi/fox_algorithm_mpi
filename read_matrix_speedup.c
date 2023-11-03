#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define dim_cart_comm 2
#define tag 0
#define VERBOSE 0
#define INF 99999

int matrixDim, *Matrix;

int dst, src, blocksquare;
MPI_Status status;
MPI_Comm gridComm, rowComm, colComm;
int myRow, myCol;


void fox_alg(int, int, int*, int*, int*, int*);

int min(int v1, int v2){
    return v1<v2 ? v1 : v2;
}
void min_plus_matrix(int* a, int* b, int* c, int blksize){
  int i,j,k;
  
  for (i = 0; i < blksize; i++)
    for (j = 0; j < blksize; j ++)
      for (k = 0; k < blksize; k++) {
        c[i*blksize+j] = min(c[i*blksize+j], (a[i*blksize+k] + b[k*blksize+j]));
      }
}

void print_matrix(int* Matrix, int dim) {
    for (int i=0; i<dim; i++) {
        for (int j=0; j<dim-1; j++) {printf("%d ", Matrix[i*dim + j]);}
        printf("%d\n", Matrix[i*dim + dim-1]);
    }
}

// Given a dimention this function reads a matrix from that dim*dim inputs
void read_matrix(int dim, int* Matrix) {
    if (VERBOSE) {printf("Reading matrix of size (%d)\n", matrixDim);}
    for (int i=0; i<dim; ++i) {
        for (int j=0; j<dim; ++j) {
            scanf("%d", &Matrix[i*dim + j]);
            if (Matrix[i*dim + j] == 0 & i!=j) {Matrix[i*dim + j] = INF;}
        }
    }
    if (VERBOSE) {printf("Finished reading input\n");}
}
void alloc_contiguous(int** V, int dim){
    (*V) = (int *)malloc(dim*dim*sizeof(int));
}
void copy_matrix(int* V1, int* V2, int dim) {
    for (int i=0; i<dim; i++){
        for (int j=0; j<dim; j++) {
            V2[i*dim + j] = V1[i*dim + j];
        }
    }
}

int main(int argc, char **argv) {
    // Create an int and a char variable
    int numprocs, rank, numblocks, blocksize;
    int *Matrix_A, *Matrix_B, *Matrix_C, *Matrix_aux;
    long startTime = time(NULL); 

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Reading matrix from INPUT only to process 0
    if (rank == 0) {
        scanf("%d", &matrixDim);
        alloc_contiguous(&Matrix, matrixDim);
        read_matrix(matrixDim, Matrix);
    }

    // Sharing initial matrix size to all processes to initialize 
    MPI_Bcast(&matrixDim, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Validating number of processes for matrix size
    if (matrixDim % (int) sqrt(numprocs) != 0 || numprocs != pow((int)sqrt(numprocs), 2)) {
        if (rank==0){
            printf("ERROR: Invalid configuration!");
            if (VERBOSE) {
                printf("ERROR: Invalid configuration!\nproc = %d | MatrixSize = %d.\nWe must have nproc=m*m and matrix_size %% m = 0.\nPlease choose one of the following nproc:\n",  numprocs, matrixDim);
                for (int i=2; i<matrixDim; i++) {
                    if (matrixDim % i == 0) {printf("%d ", i*i);}
                }
                printf("\n");
                printf("Exiting with error after %ld seconds\n", time(NULL) - startTime);
            }
            free(Matrix);
        }
        MPI_Finalize();
        return 0;
    }

    // Initialize "Matrix" in processes which are not the "master" process
    if (rank != 0) {alloc_contiguous(&Matrix, matrixDim);}


    // Sharing initial matrix size to all processes to initialize 
    MPI_Bcast(Matrix, matrixDim*matrixDim, MPI_INT, 0, MPI_COMM_WORLD);


    // Allocating memory for final result matrix
    int* MatrixResult;
    alloc_contiguous(&MatrixResult, matrixDim);


    // Calculations for array sizes
    numblocks = sqrt(numprocs);         // Find how many blocks we'll have (num_process = numblocks*numblocks)
    blocksize = matrixDim / numblocks;  // Finding the size of each block

    // Initialize M -> Variable to store all temporary blocks (A, B, C and aux)
    Matrix_A = (int *)malloc(sizeof(int) * blocksize * blocksize);
    Matrix_B = (int *)malloc(sizeof(int) * blocksize * blocksize);
    Matrix_C = (int *)malloc(sizeof(int) * blocksize * blocksize);
    Matrix_aux = (int *)malloc(sizeof(int) * blocksize * blocksize);
    if (!(Matrix_A && Matrix_B && Matrix_C && Matrix_aux)) {
        printf( "Out of memory!\n");
        free(Matrix_A);
        free(Matrix_B);
        free(Matrix_C);
        free(Matrix_aux);
        goto exit;
    }
    
    // Initializing communicator variables;
    int  sizes[dim_cart_comm] = {numblocks, numblocks};
    int   wrap[dim_cart_comm] = {1, 1};
    int colDir[dim_cart_comm] = {1, 0};
    int rowDir[dim_cart_comm] = {0, 1};
    int optm = 1, gridRank, gridCoords[dim_cart_comm], row, col;
    blocksquare = blocksize*blocksize;
    // Creating a communicator based on a topology by dividing the MPI_COMM_WORLD. Essentially, it forms a grid
    MPI_Cart_create(MPI_COMM_WORLD, dim_cart_comm, sizes, wrap, optm, &gridComm);

    // Attributing a gridRank for process inside the gridComm
    MPI_Comm_rank(gridComm, &gridRank);
    
    // "Coordinates" of the processes given their ranks.
    MPI_Cart_coords(gridComm, gridRank, dim_cart_comm, gridCoords);

    // Dividing comunicator into rows and columns
    MPI_Cart_sub(gridComm, colDir, &colComm);   // Esta a partir o comuncador em colunas
    MPI_Cart_sub(gridComm, rowDir, &rowComm);   // Esta a partir o comuncador em linhas

    // My coordinates in the Cart Grid
    myRow = gridCoords[0];
    myCol = gridCoords[1];

    src = ( myRow + 1 ) % numblocks;            // Process Source in our cart grid (neighbour)
    dst = ( myRow + numblocks - 1) % numblocks; // Process Destination in our cart grid (neighbour)

    if(VERBOSE)
      printf("Clone %d has src = %d and dst = %d\n", gridRank, src, dst);

    // Creating each process A and B block matrix
    for(int i = 0; i < blocksize; i++) {
        for (int j=0; j<blocksize; j++) {
            Matrix_A[i*blocksize + j] = Matrix[(myRow*blocksize + i)*matrixDim + myCol*blocksize + j];
            Matrix_B[i*blocksize + j] = Matrix[(myRow*blocksize + i)*matrixDim + myCol*blocksize + j];
            Matrix_C[i*blocksize + j] = INF;
        }
    }

    for (int step=1; step<matrixDim-1; step*=2) {

        // From D_k find D_2k = D_k*D_k
        fox_alg(numblocks, blocksize, Matrix_A, Matrix_B, Matrix_C, Matrix_aux); 

        // Copy Block C to Blocks A and B to repeat process
        copy_matrix(Matrix_C, Matrix_A, blocksize);
        copy_matrix(Matrix_C, Matrix_B, blocksize);
    }


    // Writing C
    MPI_Gather(Matrix_C, blocksquare, MPI_INT, MatrixResult, blocksquare, MPI_INT, 0, MPI_COMM_WORLD);

    // Switch order of matrix result
    if (rank==0) {
        int *tmp = (int *)malloc(matrixDim*matrixDim * sizeof(int));
        int next_elem = 0;
        for (int n=0; n<numblocks*numblocks; n++) {
            for (int i=0; i<blocksize; i++) {
                for (int j=0; j<blocksize; j++) {
                    row = (int) n / numblocks;
                    col = n % numblocks;
                    tmp[i*matrixDim + j + col*blocksize + row*blocksize*matrixDim] = MatrixResult[next_elem]<INF ? MatrixResult[next_elem] : 0;
                    next_elem++;
                }
            }
        }
        for (int i=0; i<matrixDim*matrixDim; i++) {
            MatrixResult[i] = tmp[i];
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) {print_matrix(MatrixResult, matrixDim);}

    // Freeing matrix memory space
    free(Matrix);
    free(MatrixResult);
    free(Matrix_A);
    free(Matrix_B);
    free(Matrix_C);
    free(Matrix_aux);
    MPI_Finalize();
    return 0;

exit:
    printf("Exiting with error after %ld seconds\n", time(NULL) - startTime);
    MPI_Finalize();
    return 0;
}

void fox_alg(int numblocks, int blocksize, int* Matrix_A, int* Matrix_B, int* Matrix_C, int* Matrix_aux) {
    int i, j, root;
    for (i = 0; i < numblocks; i++) {
        root = ( myRow  + i ) % numblocks;
        if(root == myCol){
            // snd my matrix A to neighbour
            MPI_Bcast(Matrix_A, blocksquare, MPI_INT, root, rowComm);
            min_plus_matrix(Matrix_A, Matrix_B, Matrix_C, blocksize); 
        } else {
            // rcv an A from neighbour
            MPI_Bcast(Matrix_aux, blocksquare, MPI_INT, root, rowComm);
            min_plus_matrix(Matrix_aux, Matrix_B, Matrix_C, blocksize); 
        }
        // rotate B's 
        MPI_Sendrecv_replace(Matrix_B, blocksquare, MPI_INT, dst, tag, src, tag, colComm, &status);
    }
}