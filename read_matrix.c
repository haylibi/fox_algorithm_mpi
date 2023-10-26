#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define dim_cart_comm 2
// deixar agora blocksize 2, podemos por a pedir input de um ou calcular só um fixe
#define blocksize 2
int N, *Matrix;

void mat(int r, int c, int m1[r][c], int m2[r][c], int mr[r][c]) {
    for(int i = 0; i < r; i++) {
        for(int j = 0; j < c; j++) {
            mr[i][j] = m1[i][j] + m2[i][j];
        }
    }
}

int main(int argc, char **argv) {
    // Create an int and a char variable
    int numprocs, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        scanf("%d", &N);
        fflush(stdout);
        printf("Matrix of size (%d, %d)\n", N, N);
        
        Matrix = (int *)malloc(N * N * sizeof(int));
        for (int i=0; i<N; ++i) {
            for (int j=0; j<N; ++j) {
                scanf("%d", &Matrix[N*i + j]);
            }
        }
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank != 0) {
        Matrix = (int *)malloc(N * N * sizeof(int));
    }
    
    // fox_alg(argc, argv, N), randomly placed, só para guardar

    MPI_Bcast(Matrix, N*N, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank != 0) {
        fflush(stdout);
        printf("Process %d | N=%d\n", rank, N);

        for (int i=0; i<N; ++i) {
            printf("Process %d row %d |  ", rank, i);
            for (int j=0; j<N; ++j) {
                printf("%d ", Matrix[i*N + j]);
            }
            printf("\n");
        }
    }

    free(Matrix);
    MPI_Finalize();
    exit(0);
}

int fox_alg(int argc, char** argv, int dim) {
    int numprocs, rank, matrixDim;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    matrixDim = dim * blocksize;

    MPI_Comm gridComm, rowComm, colComm; 
    int sizes[dim_cart_comm] = { dim, dim };
    int wrap[dim_cart_comm] = { 1, 1 };
    int colDir[dim_cart_comm] = { 1, 0 };
    int rowDir[dim_cart_comm] = { 0, 1 };

    // not sure para o que e o optm 1 rn, temos que ver
    int optm = 1, gridRank, gridCoords[dim_cart_comm], row, col;
    int dst, src, blocksquare = blocksize*blocksize;
     
    // Communicator creation based on a topology, we are using the comm_world and dividing it. basicamente é uma grid
    // Podemos ter que criar um novo comunicador caso não queiramos usar o comum, will see
    MPI_Cart_create(MPI_COMM_WORLD, dim_cart_comm, sizes, wrap, optm, &gridComm);
    // Processadores
    MPI_Comm_rank(gridComm, &gridRank);
    // "Coordenadas" dos processador given o rank
    MPI_Cart_coords(gridComm, gridRank, dim_cart_comm, gridCoords);

    // Está a partir o comuncador em colunas
    MPI_Cart_sub(gridComm, colDir, &colComm);

    // Está a partir o comuncador em linhas
    MPI_Cart_sub(gridComm, rowDir, &rowComm);
}
