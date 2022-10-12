ARRAY_SIZE=26 #Please note that this is not 28 elements but rather 2^28 elements
MPI_CORES=8
OMP_CORES=8
echo -------------Compiling Files-------------
export OMP_NUM_THREADS=$OMP_CORES
make
echo
echo -------Array is size: 2^$ARRAY_SIZE-------
echo
echo ------------Serial Scan is starting---------------
./scan $ARRAY_SIZE |& tee -a terminal.out
echo -------------Serial Scan is done-------------
echo
echo ------------OMP Scan is starting---------------
./scan_omp $ARRAY_SIZE |& tee -a terminal.out
echo -------------OMP Scan is done-------------
echo
echo ------------MPI Scan is starting---------------
mpirun -np $MPI_CORES ./scan_mpi $ARRAY_SIZE |& tee -a terminal.out
echo -------------MPIs Scan is done-------------
echo
echo -------------Cleaning Files-------------
make clean
rm -vf terminal.out