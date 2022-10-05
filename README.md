# Parallel Computing Assignment COMS3008A  
## 2303289, Damion Harvey 
  
  ### How to run & compile:
  Each program is setup in its own folder (ie, scan, sssp_dijsktra). Using a bash terminal, cd into the desired directory.  
  The provided shell file (run.sh) will compile all 3 programs, run all 3 files (displaying the useful info to the screen) and then discard the executables.  
  Please ensure to have the correct version of g++, OpenMP, and OpenMPI installed. 
  
  ### Examples of output for the seperate collection of programs:
  
  
### Project File Stucture is as follows: 
```
PC-Assignment
│   README.md
│   ResearchAnalysis.pdf
│
└───scan
│   │   makefile
│   │   run.sh
│   │   scan.cpp
│   │   scan_omp.cpp
│   │   scan_mpi.cpp
│   
│   
└───sssp_dijsktra
    └───graphs
    │   makefile
    │   run.sh
    │   sssp.cpp  
    │   sssp_omp.cpp      
    │   sssp_mpi.cpp       
```
  
  
  ### Other information:
  All files were written using the C++11 standard. Files were all compiled using the following versions:
  Ubuntu WSL version ubuntu1~20.04.1
  g++ --> g++ (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
  OpenMP --> (Bundled with the g++ from above)
  OpenMPI --> mpirun (Open MPI) 4.0.3
  
