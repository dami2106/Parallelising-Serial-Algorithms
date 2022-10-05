ARRAY_SIZE=28
ITERATIONS=1
echo -------------Compiling Files-------------
make
echo
echo -------Array is size: 2^$ARRAY_SIZE, with $ITERATIONS iterations-------
echo
echo ------------Serial Scan is starting---------------
SER_TIME=$(./scan $ARRAY_SIZE |& tee -a terminal.out)
echo Serial Time Taken : $SER_TIME
echo -------------Serial Scan is done-------------
echo
echo ------------OMP Scan is starting---------------
OMP_TIME=$(./scan_omp $ARRAY_SIZE |& tee -a terminal.out)
echo Speed-Up : $OMP_TIME
echo -------------OMP Scan is done-------------
echo
echo ------------MPI Scan is starting---------------
MPI_TIME=$(mpirun -np 8 ./scan_mpi $ARRAY_SIZE |& tee -a terminal.out)
echo Speed-Up : $MPI_TIME
echo -------------OMP Scan is done-------------
#echo
#echo Serial Time : $SER_TIME
#echo OMP Speed Up : $OMP_TIME
#echo MPI Speed Up : $MPI_TIME
echo
echo -------------Cleaning Files-------------
make clean
rm -vf terminal.out