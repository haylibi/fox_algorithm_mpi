#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mmlib.h"

#define dim_cart_comm 2
#define tag 0
#define VERBOSE 1
int matrixDim, **Matrix;
int* contiguousMatrix;

int fox_alg(int, int, int, int**, int, int**, int**, int**);


// Matrix multiplication function
void mat_(int r, int c, int dim, int* m1, int* m2, int* mr) {
    for(int i = 0; i < r; i++) {
        for(int j = 0; j < c; j++) {
            mr[i*dim + j] = m1[i*dim + j] + m2[i*dim + j];
        }
    }
}
void mat(int dim, int* m1, int* m2, int* mr) {
    for (int i=0; i<dim; i++) 
        for (int j=0; j<dim; j++) 
            mat_(i, j, dim, m1, m2, mr);    
}


void print_matrix(int* Matrix[], int dim) {
    for (int i=0; i<dim; i++) {
        for (int j=0; j<dim; j++) {
            fprintf(stderr, "%d ", Matrix[i][j]);
        }
        fprintf(stderr, "\n");
    }
}


// Given a dimention this function reads a matrix from that dim*dim inputs
void read_matrix(int dim, int* Matrix[]) {
    if (VERBOSE) {
        fprintf(stderr, "Reading matrix of size (%d)\n", matrixDim);
    }
    for (int i=0; i<dim; ++i) {
        for (int j=0; j<dim; ++j) {
            scanf("%d", &Matrix[i][j]);
        }
    }
    if (VERBOSE) {
        fprintf(stderr, "Finished reading input\n");
    }
}

int malloc2dfloat(float ***array, int n, int m) {

    /* allocate the n*m contiguous items */
    float *p = (float *)malloc(n*m*sizeof(float));
    if (!p) return -1;

    /* allocate the row pointers into the memory */
    (*array) = (float **)malloc(n*sizeof(float*));
    if (!(*array)) {
       free(p);
       return -1;
    }

    /* set up the pointers into the contiguous memory */
    for (int i=0; i<n; i++) 
       (*array)[i] = &(p[i*m]);

    return 0;
}

int free2dfloat(float ***array) {
    /* free the memory - the first element of the array is at the start */
    free(&((*array)[0][0]));

    /* free the pointers into the memory */
    free(*array);

    return 0;
}

int main(int argc, char **argv) {
    // Create an int and a char variable
    int numprocs, rank, numblocks, blocksize;
    int *M[4];  // A_ij = M[0]  B_hk = M[1] Cxy = M[2]  TMP  = M[3]     (These are all blocks [sub matrixes]) 
    long startTime = time(NULL); 

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Reading matrix from INPUT
    if (rank == 0) {
        scanf("%d", &matrixDim);
        if (VERBOSE) {
            fprintf(stderr, "Matrix of size (%d, %d)\n", matrixDim, matrixDim);
        }
        Matrix = (int **)malloc(matrixDim * sizeof(int *));
        for (int i=0; i<matrixDim; i++) 
            Matrix[i] = (int *)malloc(matrixDim * sizeof(int));
        read_matrix(matrixDim, Matrix);
    }

    // Sharing initial matrix size to all processes to initialize 
    MPI_Bcast(&matrixDim, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Validating number of processes for matrix size
    if (numprocs % (int) sqrt(matrixDim) != 0) {
        if (rank==0){    
            fprintf(stderr, "nproc = %d | MatrixSize = %d. We must have p=m*m and matrix_size %% m = 0.\n",  numprocs, matrixDim);
            for (int i = 0; i < matrixDim; ++i) {
                free(Matrix[i]);
            }
            free(Matrix);
            fprintf(stderr, "Exiting with error after %ld seconds\n", time(NULL) - startTime);
        }
        MPI_Finalize();
        return 0;
    }

    // Initialize "Matrix" in processes which are not the "master" process
    if (rank != 0) {
        Matrix = (int **)malloc(matrixDim * sizeof(int *));
        for (int i=0; i<matrixDim; i++) 
            Matrix[i] = (int *)malloc(matrixDim * sizeof(int));
    }

    // Sharing initial matrix size to all processes to initialize 
    for (int i=0; i<matrixDim; i++) {
        MPI_Bcast(Matrix[i], matrixDim, MPI_INT, 0, MPI_COMM_WORLD);
    }


    // Allocating memory for final result matrix
    int *MatrixResult[matrixDim];  // A_ij = M[0]  B_hk = M[1] Cxy = M[2]  TMP  = M[3]     (These are all blocks [sub matrixes]) 
    for (int i=0; i<matrixDim; i++) {
        MatrixResult[i] = (int *)malloc(matrixDim * sizeof(int));
    }

    // Calculations for array sizes
    numblocks = sqrt(numprocs);         // Find how many blocks we'll have (num_process = numblocks*numblocks)
    blocksize = matrixDim / numblocks;  // Finding the size of each block

    
    // Initializing each blocks sub matrixes (A, B, Result and a temporary one)
    for (int i=0; i<4; i++)
        M[i] = (int *)malloc(blocksize * sizeof(int));
    
    if (!(M[0] && M[1] && M[2] && M[3])) {
        fprintf(stderr,  "Out of memory!\n");
        for(int i = 0; i < 4; i++) 
            free(M[i]);
        goto exit;
    }


    if (VERBOSE) {fprintf(stderr, "Process %d: Starting Fox_alg\n", rank);}
    fox_alg(rank, numblocks, blocksize, M, matrixDim, Matrix, Matrix, MatrixResult); 
    if (VERBOSE) {fprintf(stderr, "Process %d: Finished Fox_alg\n", rank);}

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0){    
        fprintf(stderr, "Final Matrix\n");
        print_matrix(MatrixResult, matrixDim);
    }

    for (int i = 0; i < matrixDim; ++i) {free(Matrix[i]);}
    free(Matrix);
    MPI_Finalize();
    return 0;

exit:
    fprintf(stderr, "Exiting with error after %ld seconds\n", time(NULL) - startTime);
    MPI_Finalize();
    return 0;
}

int fox_alg(int rank, int numblocks, int blocksize, int* M[], int matrixDim, int** M0, int** M1, int** MR) {
    if (VERBOSE) {
        fprintf(stderr, "   <fox_alg> Process %d: Inside fox_alg. Vars: \n       NumBlocks=%d | BlockSize=%d | M[0][0]=%d | MatrixDim=%d | M0[0][0]=%d | M1[0][0]=%d | MR[0][0]=%d\n", rank, numblocks, blocksize, M[0][0], matrixDim, M0[0][0], M1[0][0], MR[0][0]);
    }
    
    int i, j, myRow, myCol;
    MPI_Comm gridComm, rowComm, colComm;

    int  sizes[dim_cart_comm] = { numblocks, numblocks };
    int   wrap[dim_cart_comm] = { 1, 1 };
    int colDir[dim_cart_comm] = { 1, 0 };
    int rowDir[dim_cart_comm] = { 0, 1 };
    int optm = 1, gridRank, gridCoords[dim_cart_comm], row, col;
    int dst, src, blocksquare = blocksize*blocksize;
    MPI_Status status;

    // Communicator creation based on a topology, we are using the comm_world and dividing it. basicamente é uma grid
    // Podemos ter que criar um novo comunicador caso não queiramos usar o comum, will see
    MPI_Cart_create(MPI_COMM_WORLD, dim_cart_comm, sizes, wrap, optm, &gridComm);

    // Attributing a gridRank for process inside the gridComm
    MPI_Comm_rank(gridComm, &gridRank);
    
    // "Coordenadas" dos processos given o rank
    MPI_Cart_coords(gridComm, gridRank, dim_cart_comm, gridCoords);

    // Dividing comunicator into rows and columns
    MPI_Cart_sub(gridComm, colDir, &colComm);   // Esta a partir o comuncador em colunas
    MPI_Cart_sub(gridComm, rowDir, &rowComm);   // Esta a partir o comuncador em linhas

    // My coordinates in the Cart Grid
    myRow = gridCoords[0];
    myCol = gridCoords[1];

    if (VERBOSE) {fprintf(stderr, "   <fox_alg> Process %d - Initial grid coords (%d,%d) | Rank %d\n", rank, myRow, myCol, gridRank);}
    
    // for(i = 0; i < blocksize; i++) {
    //     get_block_row(**M0, matrixDim, myRow, myCol, i, blocksize, &(M[0][i*blocksize]));
    //     get_block_row(**M1, matrixDim, myRow, myCol, i, blocksize, &(M[1][i*blocksize]));
    // }
    // if (VERBOSE) {
    //     fprintf(stderr, "   <fox_alg> Process %d - Matrix M[0] 1st step: ", rank);
    //     print_matrix(&M[0], blocksize);
    // }

    return 0;
}

// int fox_alg(int rank, int numblocks, int blocksize, int* M[], int matrixDim, int** M0, int** M1, int** MR) {
    
    
//     for(i = 0; i < blocksize; i++) {
//         get_block_row(**M0, matrixDim, myRow, myCol, i, blocksize, &(M[0][i*blocksize]));
//         get_block_row(**M1, matrixDim, myRow, myCol, i, blocksize, &(M[1][i*blocksize]));
//     }

    
//     src = ( myRow + 1 ) % blocksize;      // Process Source in our cart grid (neighbour)
//     dst = ( myRow + blocksize - 1) % blocksize; // Process Destinatoin in our cart grid (neighbour)

//     for (i = 0; i < blocksize; i++) {
//       int root = ( myRow  + i ) % blocksize;
//       if(root == myCol){
//         //snd my matrix A to neighbour
//         MPI_Bcast(M[0], blocksquare, MPI_INT, root, rowComm);
//         block_mult(M[2],M[0],M[1],blocksize); 
//       } else {
//         //rcv an A from neighbour
//         MPI_Bcast(M[3], blocksquare, MPI_INT, root, rowComm);
//         block_mult(M[2], M[3], M[1],blocksize); 
//       }

//       //rotate B's 
//       MPI_Sendrecv_replace(M[1], blocksquare, MPI_INT, dst, tag, src, tag, colComm, &status);

//     }

//     fprintf(stderr, "Clone %d computed:\n\n", rank);
//     for(i = 0; i < blocksquare; i++){
//       fprintf(stderr, "%5d ", (M[2][i]));
//       if((( i + 1) % blocksize) == 0)fprintf(stderr, "\n");
//     }
//     fprintf(stderr, "\n");
    
    
    
//     // Writing C
//     fprintf(stderr, "Writing C\n");
//     for(i = 0; i < blocksize; i++)
//       set_block_row(**MR, matrixDim, myRow, myCol, i, blocksize, &(M[2][i*blocksize]));

//     return 0;
// }
