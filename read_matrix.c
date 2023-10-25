#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char **argv) {
    // Create an int and a char variable
    int myNum, numprocs, rank;
    char myChar;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        printf("Please give N= ");
        fflush(stdout);
        scanf("%d", &myNum);
        // Print the number
        printf("Your number1 is: %d\n", myNum);
        fflush(stdout);
        scanf("%d", &myNum);
        printf("Your number2 is: %d\n", myNum);
        fflush(stdout);
        scanf("%d", &myNum);
        printf("Your number3 is: %d\n", myNum);
    }

    MPI_Finalize();
    exit(0);
}