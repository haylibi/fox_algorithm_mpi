# Makefile for shortest_path.c using mpicc

# MPI compiler and compiler flags
MPICC = mpicc
CFLAGS = -lm

# Source file and executable name
SRC = shortest_path.c
EXECUTABLE = fox

# Default target (executable)
$(EXECUTABLE): $(SRC)
	$(MPICC) $(SRC) -o $(EXECUTABLE)  $(CFLAGS) 

# Phony target to clean up generated files
clean:
	rm -f $(EXECUTABLE)
