# Parallel Computing Assignment COMS3008A  
## 2303289, Damion Harvey 
  
### How to run & compile:
Each program is setup in its own folder (ie, scan, sssp_dijsktra). Using a bash terminal, cd into the desired directory using the 'cd' command.  
The provided shell file (run.sh - executed wia './run.sh') will compile all 3 programs, run all 3 files (displaying the useful info to the screen) and then discard the executables.  
Please ensure to have the correct version of g++, OpenMP, and OpenMPI installed.   
*Please see examples at the bottom of the read me of running the different program sets*
  
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
  Ubuntu WSL --> version ubuntu1 20.04.1
  g++ --> g++ Ubuntu 9.4.0
  OpenMP --> (Bundled with the g++ from above)
  OpenMPI --> (Open MPI) 4.0.3
    
    
 ### Examples of output for the seperate collection of programs:
 #### Scan Example:
 <img width="613" alt="Scan Example" src="https://user-images.githubusercontent.com/102320150/194163742-772589d3-7ef0-4874-9cbc-0bdbf33d021a.png">  
   
 #### SSSP Example:
 <img width="626" alt="SSSP Example" src="https://user-images.githubusercontent.com/102320150/194164228-e268ffdd-9a7e-497a-acde-263878310c0f.png">


