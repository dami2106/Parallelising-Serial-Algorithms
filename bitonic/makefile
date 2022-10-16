P1=bitonic
P2=bitonic_omp
P3=bitonic_mpi
#INC="./inc"
#FLAGS=-I$(INC)
OMPFLAG=-fopenmp
MPICC=mpic++
CC=g++
CFLAGS = -Wall


all: $(P1) $(P2) $(P3)

$(P1): $(P1).cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(P1).cpp -o $(P1)

$(P2): $(P2).cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(P2).cpp -o $(P2)

$(P3): $(P3).cpp
	$(MPICC) $(CFLAGS)  $(P3).cpp -o $(P3)

clean:
	rm -vf $(P1) $(P2) $(P3)