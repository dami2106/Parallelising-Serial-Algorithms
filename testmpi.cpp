#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char **argv) {

    int num_processors, processor_id;
    // Initialize the MPI environment
    int fucker = 84539025;
    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
    MPI_Comm_rank(MPI_COMM_WORLD, &processor_id);
    int fuck;
    if (processor_id == 0) {
        cout << " Fucker: " <<fucker << endl;
    }
    MPI_Bcast(&fuck, 1, MPI_INT, 0, MPI_COMM_WORLD);
    cout << "Other : " << processor_id << " : " << fuck << endl;


    MPI_Finalize();

    return 0;
}