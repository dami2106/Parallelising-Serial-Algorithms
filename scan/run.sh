ARRAY_SIZE=28
ITERATIONS=1
make >/dev/null
echo -------Array is size: 2^$ARRAY_SIZE, with $ITERATIONS iterations-------
echo
echo ------------Serial Scan is starting---------------
SER_TIME=$(./scan $ARRAY_SIZE |& tee -a terminal.out)
echo -------------Serial Scan is done-------------
echo
echo ------------OMP Scan is starting---------------
OMP_TIME=$(./scan_omp $ARRAY_SIZE |& tee -a terminal.out)
echo -------------OMP Scan is done-------------
echo
echo ------------MPI Scan is starting---------------
MPI_TIME=$(mpirun -np 8 ./scan_mpi $ARRAY_SIZE |& tee -a terminal.out)
echo -------------OMP Scan is done-------------
echo
echo Serial Time : $SER_TIME
echo OMP Speed Up : $OMP_TIME
echo MPI Speed Up : $MPI_TIME
make clean >/dev/null
rm -vf terminal.out >/dev/null