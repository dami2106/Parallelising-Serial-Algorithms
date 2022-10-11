ARRAY_SIZE=28
echo -------------Compiling Files-------------
make
echo
echo -------Array is size: 2^$ARRAY_SIZE-------
echo
echo ------------Serial Scan is starting---------------
$(./scan $ARRAY_SIZE |& tee -a terminal.out)
echo -------------Serial Scan is done-------------
echo
echo ------------OMP Scan is starting---------------
$(./scan_omp $ARRAY_SIZE |& tee -a terminal.out)
echo -------------OMP Scan is done-------------
echo
echo ------------MPI Scan is starting---------------
$(mpirun -np 8 ./scan_mpi $ARRAY_SIZE |& tee -a terminal.out)
echo -------------OMP Scan is done-------------
echo
echo -------------Cleaning Files-------------
make clean
rm -vf terminal.out