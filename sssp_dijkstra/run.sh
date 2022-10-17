GRAPH=graph_7
MPI_CORES=8

echo -e "\e[1;90m━━━━━━━━━║ Compiling Files and Setup ║━━━━━━━━\e[0m"
MAKE=$(make)
echo
echo -e "\e[1;36m━━━║ Performing Dijsktra on graph $GRAPH ║━━━\e[0m"
echo
echo -e "\e[1;97m ━━━━━━━━━━║ Serial SSSP  ║━━━━━━━━━━ \e[0m"
SER=$(./sssp $GRAPH |& tee -a terminal.out)
echo -e "\e[1;92m$SER\e[0m"
echo -e "\e[1;97m ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ \e[0m"
echo
echo -e "\e[1;97m ━━━━━━━━━━║ OMP SSSP ║━━━━━━━━━━ \e[0m"
OMP=$(./sssp_omp $GRAPH |& tee -a terminal.out)
echo -e "\e[1;92m$OMP\e[0m"
echo -e "\e[1;97m ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ \e[0m"
echo
echo -e "\e[1;97m ━━━━━━━━━━║ MPI SSSP ║━━━━━━━━━━━ \e[0m"
MPI=$(mpirun -np $MPI_CORES ./sssp_mpi $GRAPH |& tee -a terminal.out)
echo -e "\e[1;92m$MPI\e[0m"
echo -e "\e[1;97m ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\e[0m"
echo
echo -e "\e[1;90m ━━━━━━━━━━║ Cleaning Files ║━━━━━━━━━━━━ \e[0m"
MAKECLEAN=$(make clean)
REM=$(rm -vf terminal.out)
echo
