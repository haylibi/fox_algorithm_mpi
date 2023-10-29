#!/bin/bash

# Assuming your MPI program is compiled as read_matrix.out
MPI_EXECUTABLE="./read_matrix.out"

# Input file
INPUT_DIR="./inputs"

# Output file
OUTPUT_DIR="./outputs"

start_time=$(date +"%s")

# Run the MPI program using mpirun
echo "mpirun -n 4 \"$MPI_EXECUTABLE\" < \"$INPUT_DIR/input6\" > \"$OUTPUT_DIR/output6.txt\" 2>&1"
mpirun -n 4 "$MPI_EXECUTABLE" < "$INPUT_DIR/input6" > "$OUTPUT_DIR/output6.txt" 2>&1
echo "Elapsed Time: $(($(date +"%s") - start_time)) seconds"

start_time=$(date +"%s")
echo "mpirun -n 25 \"$MPI_EXECUTABLE\" < \"$INPUT_DIR/input300\" > \"$OUTPUT_DIR/output300.txt\" 2>&1"
mpirun -n 25 "$MPI_EXECUTABLE" < "$INPUT_DIR/input300" > "$OUTPUT_DIR/output300.txt" 2>&1
echo "Elapsed Time: $(($(date +"%s") - start_time)) seconds"

# Print a message indicating the completion of the run
echo "Completed run. Output saved to $OUTPUT_DIR"
