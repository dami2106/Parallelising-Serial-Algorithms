ARRAY_SIZE=28 #Please note that this is not 28 elements but rather 2^28 elements
MPI_CORES=8

echo -e "\e[1;90m━━━━━━━━━║ Compiling Files and Setup ║━━━━━━━━\e[0m"
MAKE=$(make)
echo
echo -e "\e[1;36m━━━━━━━━━║ Array is size: 2^$ARRAY_SIZE ║━━━━━━━━━\e[0m"
echo
echo -e "\e[1;97m━━━━━━━━━━║ Serial Scan  ║━━━━━━━━━━ \e[0m"
SER=$(./scan $ARRAY_SIZE |& tee -a terminal.out)
echo -e "\e[1;92m$SER\e[0m"
echo -e "\e[1;97m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ \e[0m"
echo
echo -e "\e[1;97m━━━━━━━━━━║ OMP Scan  ║━━━━━━━━━━ \e[0m"
OMP=$(./scan_omp $ARRAY_SIZE |& tee -a terminal.out)
echo -e "\e[1;92m$OMP\e[0m"
echo -e "\e[1;97m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ \e[0m"
echo
echo -e "\e[1;97m━━━━━━━━━━║ MPI Scan  ║━━━━━━━━━━ \e[0m"
MPI=$(mpirun -np $MPI_CORES ./scan_mpi $ARRAY_SIZE |& tee -a terminal.out)
echo -e "\e[1;92m$MPI\e[0m"
echo -e "\e[1;97m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ \e[0m"
echo
echo -e "\e[1;90m━━━━━━━━━━║ Cleaning Files ║━━━━━━━━━━━━ \e[0m"
MAKECLEAN=$(make clean)
REM=$(rm -vf terminal.out)
echo