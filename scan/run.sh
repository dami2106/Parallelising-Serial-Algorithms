ARRAY_SIZE=10
ITERATIONS=5
echo -------Array is size: $ARRAY_SIZE, with $ITERATIONS iterations-------
echo
echo ------------Serial Scan is starting---------------
SER_TIME=$(./scan $ARRAY_SIZE $ITERATIONS |& tee -a terminal.out)
echo -------------Serial Scan is done-------------
echo
echo ------------OMP Scan is starting---------------
PAR_TIME=$(./scan_omp $ARRAY_SIZE $ITERATIONS |& tee -a terminal.out)
echo -------------OMP Scan is done-------------
