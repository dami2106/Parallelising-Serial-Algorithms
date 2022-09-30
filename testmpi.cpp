#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {

    int num_processors, processor_id;
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
    MPI_Comm_rank(MPI_COMM_WORLD, &processor_id);

    std::cout << "Hi from Process " << processor_id << " out of " << num_processors << std::endl;

    MPI_Finalize();

    return 0;
}