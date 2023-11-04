#!/bin/bash
#!/bin/bash

MPI_SCRIPT="read_matrix_speedup.c"
MPI_EXECUTABLE="./a.out"

# Input file
INPUT_DIR="./inputs"
# Function to get the current time in milliseconds
current_time_in_ms() {
  echo $(($(date +%s%N) / 1000000))
}


# Output file
OUTPUT_DIR="./outputs"

# Compile the C program
mpicc "$MPI_SCRIPT" -o "$MPI_EXECUTABLE" -lm


# Number of process to iterate through
numprocs=("1" "4" "9" "16" "25" "36")

# Iterate through the array
for file in "$INPUT_DIR"/*
do
    echo "############################################################"
    echo ""
    echo "                       FILE '$file'                         "
    echo ""
    echo "############################################################"
    for i in "${numprocs[@]}"
    do
        echo "mpirun --hostfile hostfile -n $i \"$MPI_EXECUTABLE\" < \"$file\" > \"$OUTPUT_DIR/output_$(basename $file)-np_$i.txt\" 2>&1"
        start_time=$(current_time_in_ms)
        mpirun --hostfile hostfile -n $i "$MPI_EXECUTABLE" < "$file" > "$OUTPUT_DIR/output_$(basename $file)-np_$i.txt" 2>&1
        echo "Elapsed Time: $((($(current_time_in_ms) - start_time))) milliseconds"
    done
done
