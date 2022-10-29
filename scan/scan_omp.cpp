/*
 * Parallel implementation of Blelloch Parallel Scan
 */
#include <iostream>
#include <vector>
#include <random>
#include <omp.h>
#include <array>
#include <algorithm>

#define NUMTHREADS 8

std::vector<int> generateArray(int N);

void serialFullScan(std::vector<int> &in, std::vector<int> &out, int N);

void ompFullScan(std::vector<int> &in, int N);

void checkInput(std::string arrSize);

void ompBlelloch(std::vector<int> &in, int N);


int main(int argc, char *argv[]) {

    //Check the input array size is correct
    checkInput(argv[1]);

    //Get the size of the array from the parameters
    int N = (int) pow(2, atoi(argv[1]));
    //Set the variables used for timing
    double startTime, sRunTime = 0, pRunTime = 0, bleRunTime = 0;

    //Generate an array of random elements of size N and back it up
    std::vector<int> in = generateArray(N);
    std::vector<int> blIn = in;

    //An array of 0s of size N used to store the result of serial full scan
    std::vector<int> out(N, 0);

    //Time and call the serial full scan
    startTime = omp_get_wtime();
    serialFullScan(in, out, N);
    sRunTime = omp_get_wtime() - startTime;

    //Time and call the classic blelloch
    startTime = omp_get_wtime();
    ompBlelloch(in, N);
    pRunTime = omp_get_wtime() - startTime;

    //Time and call the naive blelloch
    startTime = omp_get_wtime();
    ompFullScan(blIn, N);
    bleRunTime = omp_get_wtime() - startTime;

    //Validate the parallel data against the serial sum array and serial parallel array
    if ((in != out) || (blIn != out)) {
        std::cout << "(Validation Failed!)\n";
    } else {
        std::cout << "(Validation Passed!)\n";
        std::cout << "Serial Time : " << sRunTime << std::endl << std::endl;
        std::cout << "Classic Blelloch Parallel Time : " << pRunTime << std::endl;
        std::cout << "Classic Blelloch Speed-Up : " << (sRunTime / pRunTime) << std::endl << std::endl;
        std::cout << "Naive Blocking Blelloch Parallel Time : " << bleRunTime << std::endl;
        std::cout << "Naive Blocking Blelloch Speed-Up : " << (sRunTime / bleRunTime) << std::endl;
    }

    return 0;
}

/*
 * Checks if the provided argument is correct, exits if it isn't
 */
void checkInput(std::string arrSize) {
    if (pow(2, stoi(arrSize)) < NUMTHREADS) std::cout << "INCORRECT ARRAY SIZE\n", exit(1);
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
 * A function that performs a work efficient blelloch scan on the array
 */
void ompBlelloch(std::vector<int> &in, int N) {
    //Initialise variables needed in the algorithm
    int d, k, inc, ind1, ind2, i, t, temp, inc2;
    //Create the parallel region
#pragma omp parallel num_threads(NUMTHREADS) shared(in, temp) private(d, k, inc, ind1, ind2, i, t, inc2)
    {
        //Perform the up sweep(reduction) on the data by referencing a tree in memory
        //And adding up paired elements
        //Here the tree is traversed from the root node to the leaf nodes
        for (d = 0; d < (int) log2(N); d++) {
            inc = (1 << (d + 1)); //Bit shift inc , same as writing 2^(d+1) just faster
            inc2 = (1 << d); //Bit shift inc2, same as writing 2^d, just faster
#pragma omp for
            for (k = 0; k < N - 1; k += inc) {
                ind1 = k + inc - 1;
                ind2 = k + inc2 - 1;
                in[ind1] += in[ind2];
            }
        }
        //The root node of the tree now holds the max sum value of the array

        //We then back up the end element and set it to 0 in order to do a pre scan
#pragma omp single
        {
            temp = in[N - 1];
            in[N - 1] = 0;
        }


//We now perform down sweep on the tree,
        for (d = log2(N) - 1; d >= 0; --d) {
            inc = (1 << (d + 1));
            inc2 = (1 << d);
#pragma omp for
            for (i = 0; i <= N - 1; i += inc) {
                t = in[i + inc2 - 1];
                in[i + inc2 - 1] = in[i + inc - 1];
                in[i + inc - 1] = t + in[i + inc - 1];
            }
        }

        //We then need to convert the pre-scan to the inclusive scan by "shifting" to the left and restoring the max
#pragma omp single
        {
            in.erase(in.begin());
            in.emplace_back(temp);
        }
    }
}


/*
 * A function that performs a naive blelloch scan on the array
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