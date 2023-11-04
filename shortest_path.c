#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define dim_cart_comm 2
#define tag 0
#define VERBOSE 0
#define print_execution_time 0
#define INF 999999

int matrixDim, *Matrix;

int dst, src, blocksquare;
MPI_Comm gridComm, rowComm, colComm;
int myRow, myCol;

// Initializing Fox algorithm
void fox_alg(int, int, int*, int*, int*, int*);

// Defining Min-Plus calculation between two Matrices
int min(int v1, int v2){return v1<v2 ? v1 : v2;}
void min_plus_matrix(int* a, int* b, int* c, int blksize){
  int i,j,k;
  
  for (i = 0; i < blksize; i++)
    for (j = 0; j < blksize; j ++)
      for (k = 0; k < blksize; k++) {
        c[i*blksize+j] = min(c[i*blksize+j], (a[i*blksize+k] + b[k*blksize+j]));
      }
}

// Function which prints out a contiguous Matrix
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
            // Replacing 0's with "INF" in order to calculate correctly the "Min-Plus"
            if ((Matrix[i*dim + j] == 0) & (i!=j)) {Matrix[i*dim + j] = INF;}
        }
    }
    if (VERBOSE) {printf("Finished reading input\n");}
}

// Allocates contiguous memory for a 2D array of integers with dimensions dim x dim
void alloc_contiguous(int** V, int dim){
    (*V) = (int *)malloc(dim*dim*sizeof(int));
}

// Coppies matrix V1 to V2 and both arrays are assumed to have sizes dim x dim
void copy_matrix(int* V1, int* V2, int dim) {
    for (int i=0; i<dim; i++){
        for (int j=0; j<dim; j++) {
            V2[i*dim + j] = V1[i*dim + j];
        }
    }
}


int main(int argc, char **argv) {
    int numprocs, rank, numblocks, blocksize;
    int *Matrix_A, *Matrix_B, *Matrix_C, *Matrix_aux;
    double start;

    // Initializing MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    // Reading matrix from INPUT only to process 0
    if (rank == 0) {
        scanf("%d", &matrixDim);
        alloc_contiguous(&Matrix, matrixDim);
        read_matrix(matrixDim, Matrix);
    }
    start = MPI_Wtime(); // Start counting time after reading input matrix

    // Sharing initial matrix size to all processes to initialize 
    MPI_Bcast(&matrixDim, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Validating number of processes for matrix size
    if (matrixDim % (int) sqrt(numprocs) != 0 || numprocs != pow((int)sqrt(numprocs), 2)) {
        if (rank==0){
            printf("ERROR: Invalid configuration!");
            free(Matrix);
        }
        MPI_Finalize();
        return 0;
    }

    // Initialize "Matrix" in processes which are not the "master" process
    if (rank != 0) {alloc_contiguous(&Matrix, matrixDim);}
    // Allocating memory for final result matrix
    int* MatrixResult;
    alloc_contiguous(&MatrixResult, matrixDim);


    // Sharing initial matrix size to all processes to initialize 
    MPI_Bcast(Matrix, matrixDim*matrixDim, MPI_INT, 0, MPI_COMM_WORLD);


    // Calculations for array sizes
    numblocks = sqrt(numprocs);         // Find how many blocks we'll have (num_process = numblocks*numblocks)
    blocksize = matrixDim / numblocks;  // Finding the size of each block
    // Initializing communicator variables;
    int  sizes[dim_cart_comm] = {numblocks, numblocks};
    int   wrap[dim_cart_comm] = {1, 1};
    int colDir[dim_cart_comm] = {1, 0};
    int rowDir[dim_cart_comm] = {0, 1};
    int optm = 1, gridRank, gridCoords[dim_cart_comm], row, col;
    blocksquare = blocksize*blocksize;
    // Communicator creation based on a topology, we are using the comm_world and dividing it. basicamente é uma grid
    // Podemos ter que criar um novo comunicador caso não queiramos usar o comum, will see
    MPI_Cart_create(MPI_COMM_WORLD, dim_cart_comm, sizes, wrap, optm, &gridComm);
    MPI_Comm_rank(gridComm, &gridRank);     // Attributing a gridRank for process inside the gridComm
    MPI_Cart_coords(gridComm, gridRank, dim_cart_comm, gridCoords);  // "Coordenadas" dos processos given o rank
    MPI_Cart_sub(gridComm, colDir, &colComm);   // Esta a partir o comuncador em colunas
    MPI_Cart_sub(gridComm, rowDir, &rowComm);   // Esta a partir o comuncador em linhas

    // Process coordinates in the Cart Grid
    myRow = gridCoords[0];
    myCol = gridCoords[1];

    src = ( myRow + 1 ) % numblocks;            // Process Source in our cart grid (neighbour)
    dst = ( myRow + numblocks - 1) % numblocks; // Process Destination in our cart grid (neighbour)

    if(VERBOSE) {
      printf("Clone %d has src = %d and dst = %d\n", gridRank, src, dst);
    }


    // Initialize Matrixes for processes calculatoins
    Matrix_A   = (int *)malloc(sizeof(int) * blocksize * blocksize);
    Matrix_B   = (int *)malloc(sizeof(int) * blocksize * blocksize);
    Matrix_C   = (int *)malloc(sizeof(int) * blocksize * blocksize);
    Matrix_aux = (int *)malloc(sizeof(int) * blocksize * blocksize);
    if (!(Matrix_A && Matrix_B && Matrix_C && Matrix_aux)) {
        printf( "Out of memory!\n");
        free(Matrix_A);
        free(Matrix_B);
        free(Matrix_C);
        free(Matrix_aux);
        goto exit;
    }
    
    // Initializing matrices with values from Adjacency matrix
    for(int i = 0; i < blocksize; i++) {
        for (int j=0; j<blocksize; j++) {
            Matrix_A[i*blocksize + j] = Matrix[(myRow*blocksize + i)*matrixDim + myCol*blocksize + j];
            Matrix_B[i*blocksize + j] = Matrix[(myRow*blocksize + i)*matrixDim + myCol*blocksize + j];
            Matrix_C[i*blocksize + j] = Matrix[(myRow*blocksize + i)*matrixDim + myCol*blocksize + j];
        }
    }

    // Calculation of Shortest Paths
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
    // (Since gather returns a contiguous matrix and joins processes from left to right up to down, we change that to correct order)
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
    

    // Only print Execution time if requested
    if (print_execution_time & (rank==0)) {
        printf("Execution time: %.5lf\n", MPI_Wtime() - start);
    }

    // Print matrix out
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
    printf("Exiting with error after %.5lf seconds\n", MPI_Wtime() - start);
    MPI_Finalize();
    return 0;
}

void fox_alg(int numblocks, int blocksize, int* Matrix_A, int* Matrix_B, int* Matrix_C, int* Matrix_aux) {
    int i, root;
    // Iterating through the number of blocks
    for (i = 0; i < numblocks; i++) {
        root = ( myRow  + i ) % numblocks; // Finding root matrix in block (root starts on the diagonals)
        if(root == myCol){
            //snd my matrix A to neighbour
            MPI_Bcast(Matrix_A, blocksquare, MPI_INT, root, rowComm);
            min_plus_matrix(Matrix_A, Matrix_B, Matrix_C, blocksize); 
        } else {
            //rcv an A from neighbour
            MPI_Bcast(Matrix_aux, blocksquare, MPI_INT, root, rowComm);
            min_plus_matrix(Matrix_aux, Matrix_B, Matrix_C, blocksize); 
        }
        //   rotate B's 
        MPI_Sendrecv_replace(Matrix_B, blocksquare, MPI_INT, dst, tag, src, tag, colComm, MPI_STATUS_IGNORE);
    }
}
