P1=scan
P2=scan_omp

#INC="./inc"
#FLAGS=-I$(INC)
OMPFLAG=-fopenmp
MPIFLAG=-np 16
MPICC=mpic++
CC=g++
CFLAGS = -Wall


all: $(P1) $(P2)

$(P1): $(P1).cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(P1).cpp -o $(P1)

$(P2): $(P2).cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(P2).cpp -o $(P2)

clean:
	rm -vf $(P1) $(P2)