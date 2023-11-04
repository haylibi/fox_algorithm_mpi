# fox_algorithm_mpi
Implementation of the Fox algorithm using OpenMPI Parallel Programming

Report Overleaf URL:
[https://www.overleaf.com/7311244184fhkjshjzvdpy#c0d1fc](https://www.overleaf.com/7311244184fhkjshjzvdpy#c0d1fc)

# How to run
1. Compile the program "read_matrix_speedup.c" with the following mpi command `mpicc read_matrix_speedrun.c -lm -o a.out`
2. Running the program:
  1. Create a file named "hotsfile" containing the following content: `localhost slots=<max_num_processes_allowed>`
  2. Execute the program using the following command `mpirun --hostfile hostfile -n #numprocesses# a.out < "inputs/input6"`
