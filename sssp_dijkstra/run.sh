GRAPH=graph_7
echo -------------Compiling Files-------------
make
echo
echo -------Performing Dijsktra on graph $GRAPH-------
echo
echo ------------Serial SSSP is starting---------------
SER_TIME=$(./sssp $GRAPH |& tee -a terminal.out)
echo Serial Time Taken : $SER_TIME
echo -------------Serial SSSP is done-------------
echo
echo ------------OMP SSSP is starting---------------
OMP_TIME=$(./sssp_omp $GRAPH |& tee -a terminal.out)
echo Speed-Up : $OMP_TIME
echo -------------OMP SSSP is done-------------
echo
echo ------------MPI SSSP is starting---------------
MPI_TIME=$(mpirun -np 8 ./sssp_mpi $GRAPH |& tee -a terminal.out)
echo Speed-Up : $MPI_TIME
echo -------------OMP SSSP is done-------------
#echo
#echo Serial Time : $SER_TIME
#echo OMP Speed Up : $OMP_TIME
#echo MPI Speed Up : $MPI_TIME
echo
echo -------------Cleaning Files-------------
make clean
rm -vf terminal.out