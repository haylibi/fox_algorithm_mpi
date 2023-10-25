#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

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