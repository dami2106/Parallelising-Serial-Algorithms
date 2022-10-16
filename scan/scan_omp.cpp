/*
 * Parallel implementation of Blelloch Parallel Scan
 */
#include <iostream>
#include <vector>
#include <random>
#include <omp.h>
#include <array>

#define NUMTHREADS 8

std::vector<int> generateArray(int N);
void serialFullScan(std::vector<int> &in, std::vector<int> &out, int N);
void ompFullScan(std::vector<int> &in, int N);

int main(int argc, char *argv[]) {
    //Get the size of the array from the parameters
    int N = (int) pow(2, atoi(argv[1]));
    //Set the variables used for timing
    double startTime, sRunTime = 0, pRunTime = 0;

    //Generate an array of random elements of size N
    std::vector<int> in = generateArray(N);
    //An array of 0s of size N used to store the result of serial full scan
    std::vector<int> out(N, 0);

    //Time and call the serial full scan
    startTime = omp_get_wtime();
    serialFullScan(in, out, N);
    sRunTime += omp_get_wtime() - startTime;

    //Time and call the parallel full scan
    startTime = omp_get_wtime();
    ompFullScan(in, N);
    pRunTime += omp_get_wtime() - startTime;

    //Validate the parallel data against the serial sum array and serial parallel array
    if (in != out) {
        std::cout << "(Validation Failed!)\n";
    } else {
        std::cout << "(Validation Passed!)\n";
        std::cout << "Serial Time : " << sRunTime << std::endl;
        std::cout << "Parallel Time : " << pRunTime << std::endl;
        std::cout << "Speed-Up : " << (sRunTime / pRunTime) << std::endl;
        std::cout << "Efficiency : " << (sRunTime / pRunTime)/NUMTHREADS << std::endl;
    }

    return 0;
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
    std::uniform_int_distribution<> dist(1, 50);

    //Insert the random nubmers into the array
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
 * A function that performs a parallel full scan on the given array
 */
void ompFullScan(std::vector<int> &in, int N) {
    //Initialise variables needed for the parallel execution
    int threadID, threadCount, threadBoundLeft, threadBoundRight, i;
    //Initialise 2 arrays needed to propagate the local sums to other threads
    std::vector<int> globalSum;
    std::vector<int> incrementValues;

    //Start the parallel region
#pragma omp parallel num_threads(NUMTHREADS) private(i, threadCount, threadID, threadBoundLeft, threadBoundRight) shared(in, N, globalSum, incrementValues)
    {
        threadID = omp_get_thread_num(); //Stores the current threads rank
        threadCount = omp_get_num_threads(); //Stores the thread count

        //Sets the left & right bound of the current thread (used to ensure each thread only operates on its portion of the
        //distance array and adjacency list). Same effect as creating local arrays without the overhead (1D blocking)
        threadBoundLeft = threadID * (N / threadCount);
        threadBoundRight = ((threadID + 1) * (N / threadCount)) - 1;

//Ensure that the vectors can be indexed without initialising values
#pragma omp single
        {
            globalSum.resize(threadCount);
            incrementValues.resize(threadCount);
        }
//Synchronise the threads
#pragma omp barrier

//Compute the full scan of the sub array of each thread using the normal serial method
        for (i = threadBoundLeft + 1; i <= threadBoundRight && i < N; i++)
            in[i] += in[i - 1];
        globalSum[threadID] = in[i - 1];

//Synchronise the threads
#pragma omp barrier
//Set the sum array that stores the cumulative sum of each sub array
//Use a shift operator on i ( same as incrementing by 2^d, but faster )
        for (i = 1; i < threadCount; i <<= 1) {
            if (threadID >= i) {
                incrementValues[threadID] = globalSum[threadID] + globalSum[threadID - i];
            }
//Synchronise the threads
#pragma omp barrier
//Copy the elements to the global sums array excluding the first element
#pragma omp single
            std::copy(std::begin(incrementValues) + 1, std::end(incrementValues), std::begin(globalSum) + 1);
        }
//Synchronise the threads
#pragma omp barrier
//Add the cumulative sum to each local thread while also inserting the local values back into the initial array
        for (i = threadBoundLeft; i <= threadBoundRight; i++) {
            in[i] += globalSum[threadID] - in[threadBoundRight];
        }
    }

}