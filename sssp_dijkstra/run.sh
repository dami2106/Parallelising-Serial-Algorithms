GRAPH=graph_7
MPI_CORES=8
OMP_CORES=8
echo -------------Compiling Files and Setup-----------
export OMP_NUM_THREADS=$OMP_CORES
make
echo -----------------------------------------
echo
echo -------Performing Dijsktra on graph $GRAPH-------
echo
echo ------------Serial SSSP is starting---------------
./sssp $GRAPH |& tee -a terminal.out
echo -------------Serial SSSP is done-------------
echo
echo ------------OMP SSSP is starting---------------
./sssp_omp $GRAPH |& tee -a terminal.out
echo -------------OMP SSSP is done-------------
echo
echo ------------MPI SSSP is starting---------------
mpirun -np $MPI_CORES ./sssp_mpi $GRAPH |& tee -a terminal.out
echo -------------OMP SSSP is done-------------
echo
echo -------------Cleaning Files-------------
make clean
rm -vf terminal.out
echo -----------------------------------------
echo