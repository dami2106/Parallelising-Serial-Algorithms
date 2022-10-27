/*
 * Parallel distributed implementation of Blelloch Parallel Scan
 */
#include <iostream>
#include <vector>
#include <random>
#include <mpi.h>

std::vector<int> generateArray(int N);
void serialFullScan(std::vector<int> &in, std::vector<int> &out, int N);
void mpiFullScan(std::vector<int> &in, int N);
void checkInput(std::string arrSize, int threadCount);

int main(int argc, char *argv[]) {
    //Get the size of the array from the parameters
    int N = (int) pow(2, atoi(argv[1]));
    int id, P; //Stores the ID of the current thread

    //Set the variables used for timing
    double startTime, serRuntime = 0, parRuntime = 0;

    //Set up an array to store the initial random array
    std::vector<int> in;
    //Array to store the serial output
    std::vector<int> ser(N, 0);

    //Start MPI
    MPI_Init(NULL, NULL);

    //Get the current threads rank
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    //Only need 1 thread at a time to do the setup
    if (id == 0) {
        //Check the input array size is correct
        checkInput(argv[1], P);

        //Generate the randomised array
        in = generateArray(N);

        //Time and do the parallel fullScan
        startTime = MPI_Wtime();
        serialFullScan(in, ser, N);
        serRuntime = MPI_Wtime() - startTime;

    }

    //Synchronise all the threads
    MPI_Barrier(MPI_COMM_WORLD);

    //Start the parallel timer using 1 thread
    if(id == 0)
        startTime = MPI_Wtime();

    //Call the parallel full scan implementation
    mpiFullScan(in, N);

    //Stop the timer
    if(id == 0)
        parRuntime = MPI_Wtime() - startTime;

    if (id == 0) {
        //Validate the parallel data against the serial sum array and serial parallel array
        if (in != ser) {
            std::cout << "(Validation Failed!)\n";
        } else {
            std::cout << "(Validation Passed!)\n";
            std::cout << "Serial Time : " << serRuntime << std::endl;
            std::cout << "Parallel Time : " << parRuntime << std::endl;
            std::cout << "Speed-Up : " << (serRuntime / parRuntime) << std::endl;
            std::cout << "Efficiency : " << (serRuntime / parRuntime)/P << std::endl;
        }
    }

    //Finalise the MPI call
    MPI_Finalize();

    return 0;
}

/*
 * Checks if the provided argument is correct, exits if it isn't
 */
void checkInput(std::string arrSize, int threadCount) {
    if (pow(2, stoi(arrSize)) < threadCount) std::cout << "INCORRECT ARRAY SIZE\n", exit(1);
    if (stoi(arrSize) < 3 || stoi(arrSize) > 28) std::cout << "INCORRECT ARRAY SIZE\n", exit(1);
}

/*
 * Generate an array of size N with random elements between [1, 50)
 */
std::vector<int> generateArray(int N) {
    std::vector<int> arr(N, 0);
    //Create a blank vector of size N
    //Initialise a random device to randomly generate numbers to insert into the array
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<> dist(0, 10);

    //Insert the random numbers into the array
    for (int i = 0; i < N; i++)
        arr[i] = dist(generator);

    return arr;
}

/*
 * A function that performs a serial full scan on the given array
 */
void serialFullScan(std::vector<int> &in, std::vector<int> &out, int N) {
    out[0] = in[0]; //Set the first elements to be equal
    //Perform serial full scan
    for (int i = 1; i < N; i++)
        out[i] = in[i] + out[i - 1];
}

/*
 * A function that performs a distributed parallel full scan on the given array
 */
void mpiFullScan(std::vector<int> &in, int N) {
    //Initialise variables for each thread that will be used in the algorithm
    int threadCount, threadID, localN, localSum = 0, localIncrement = 0;

    //Get the current thread rank and total thread count
    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCount);

    //Get the size of the local array (size of the array each thread will have)
    localN = N / threadCount;

    //Initialise 2 arrays that will be used to store the full scan sums for each thread
    std::vector<int> globalSum(threadCount, 0);
    std::vector<int> localIn(localN, 0);

    //Scatter / distribute the main array to each thread
    MPI_Scatter(in.data(), localN, MPI_INT, localIn.data(), localN, MPI_INT, 0, MPI_COMM_WORLD);
//
//    //Compute the full scan sum for each thread on its subset of elements
//    for (int i = 1; i < localN; i++) {
//        localIn[i] += localIn[i - 1];
//    }
    int inc, t;
    for (int d = 0; d < (int) log2(localN); d++) {
        inc = (int) pow(2, d + 1);
        for (int k = 0; k < localN - 1; k += inc) {
            int ind1 = k + inc - 1;
            int ind2 = k + (int) pow(2, d) - 1;
            if (ind1 < localN && ind2 < localN)
                localIn[ind1] += localIn[ind2];

        }
    }


    int temp = localIn[localN - 1];
    localIn[localN - 1] = 0;

    for (int d = log2(localN) - 1; d >= 0; --d) {
        inc = (int) pow(2, d + 1);
        for (int i = 0; i <= localN - 1; i += inc) {
            t = localIn[i + pow(2, d) - 1];
            localIn[i + pow(2, d) - 1] = localIn[i + pow(2, d + 1) - 1];
            localIn[i + pow(2, d + 1) - 1] = t + localIn[i + pow(2, d + 1) - 1];
        }
    }

    localIn.push_back(temp);
    localIn.erase(localIn.begin());

    //Set the cumulative end sum for each thread
    localSum = localIn[localN - 1];

    //Gather the local sums (end values of the full scan) into the globalSums array to all threads
    MPI_Allgather(&localSum, 1, MPI_INT, globalSum.data(), 1, MPI_INT, MPI_COMM_WORLD);

    //Need to sum up the globalSums for each thread's position in the main vector
    //Ensures the local array at position i gets the sum from the end element of array (i-1)
    for (int i = 0; i < threadID; i++) {
        localIncrement += globalSum[i];
    }

    //For all local sums after the first one, add on the increment value
    if (localIncrement > 0) {
        for (int i = 0; i < localN; i++)
            localIn[i] += localIncrement;
    }

    //Finally gather all the local array sums back into the original array for thread 0
    MPI_Gather(localIn.data(), localN, MPI_INT, in.data(), localN, MPI_INT, 0, MPI_COMM_WORLD);
}
