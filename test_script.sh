#!/bin/bash

# Assuming your MPI program is compiled as read_matrix.out
MPI_SCRIPT="read_matrix_speedup.c"
MPI_EXECUTABLE="./a.out"

# Input file
INPUT_DIR="./inputs"

# Output file
OUTPUT_DIR="./outputs"

start_time=$(date +"%s")

# Compile the C program
mpicc "$MPI_SCRIPT" -o "$MPI_EXECUTABLE" -lm

# Run the MPI program using mpirun
echo "mpirun -n 4 \"$MPI_EXECUTABLE\" < \"$INPUT_DIR/input6\" > \"$OUTPUT_DIR/output6.txt\" 2>&1"
mpirun -n 4 "$MPI_EXECUTABLE" < "$INPUT_DIR/input6" > "$OUTPUT_DIR/output6.txt" 2>&1
echo "Elapsed Time: $(($(date +"%s") - start_time)) seconds"

start_time=$(date +"%s")
echo "mpirun -n 4 \"$MPI_EXECUTABLE\" < \"$INPUT_DIR/input300\" > \"$OUTPUT_DIR/output300.txt\" 2>&1"
mpirun -n 4 "$MPI_EXECUTABLE" < "$INPUT_DIR/input300" > "$OUTPUT_DIR/output300.txt" 2>&1
echo "Elapsed Time: $(($(date +"%s") - start_time)) seconds"

start_time=$(date +"%s")
echo "mpirun -n 25 \"$MPI_EXECUTABLE\" < \"$INPUT_DIR/input600\" > \"$OUTPUT_DIR/output600.txt\" 2>&1"
mpirun -n 25 "$MPI_EXECUTABLE" < "$INPUT_DIR/input600" > "$OUTPUT_DIR/output600.txt" 2>&1
echo "Elapsed Time: $(($(date +"%s") - start_time)) seconds"

start_time=$(date +"%s")
echo "mpirun -n 25 \"$MPI_EXECUTABLE\" < \"$INPUT_DIR/input900\" > \"$OUTPUT_DIR/output900.txt\" 2>&1"
mpirun -n 25 "$MPI_EXECUTABLE" < "$INPUT_DIR/input900" > "$OUTPUT_DIR/output900.txt" 2>&1
echo "Elapsed Time: $(($(date +"%s") - start_time)) seconds"

start_time=$(date +"%s")
echo "mpirun -n 36 \"$MPI_EXECUTABLE\" < \"$INPUT_DIR/input1200\" > \"$OUTPUT_DIR/output1200.txt\" 2>&1"
mpirun -n 36 "$MPI_EXECUTABLE" < "$INPUT_DIR/input1200" > "$OUTPUT_DIR/output1200.txt" 2>&1
echo "Elapsed Time: $(($(date +"%s") - start_time)) seconds"

# Print a message indicating the completion of the run
echo "Completed run. Output saved to $OUTPUT_DIR"
