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

void checkInput(const std::string &arrSize, int threadCount);

void mpiBlellochScan(std::vector<int> &in, int N);

int main(int argc, char *argv[]) {
    //Get the size of the array from the parameters
    int N = (int) pow(2, atoi(argv[1]));
    int id, P; //Stores the ID of the current thread and number of threads

    //Set the variables used for timing
    double startTime, serRuntime = 0, parRuntime = 0, bleRunTime = 0;

    //Set up an array to store the initial random arrays
    std::vector<int> in;
    std::vector<int> blellochIn;

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
        blellochIn = in;

        //Time and do the parallel fullScan
        startTime = MPI_Wtime();
        serialFullScan(in, ser, N);
        serRuntime = MPI_Wtime() - startTime;

    }

    //Synchronise all the threads
    MPI_Barrier(MPI_COMM_WORLD);

    //Start the parallel timer using 1 thread
    if (id == 0)
        startTime = MPI_Wtime();

    //Call the parallel full scan implementation
    mpiFullScan(in, N);

    //Stop the timer
    if (id == 0)
        parRuntime = MPI_Wtime() - startTime;

    //Start the parallel timer using 1 thread
    if (id == 0)
        startTime = MPI_Wtime();

    //Call the parallel full scan implementation
    mpiBlellochScan(blellochIn, N);

    //Stop the timer
    if (id == 0)
        bleRunTime = MPI_Wtime() - startTime;


    if (id == 0) {
        //Validate the parallel data against the serial sum array and serial parallel array
        if ((in != ser) || (blellochIn != ser)) {
            std::cout << "(Validation Failed!)\n";
        } else {
            std::cout << "(Validation Passed!)\n";
            std::cout << "Serial Time : " << serRuntime << std::endl << std::endl;
            std::cout << "Classic Blelloch Parallel Time : " << bleRunTime << std::endl;
            std::cout << "Classic Blelloch Speed-Up : " << (serRuntime / bleRunTime) << std::endl << std::endl;
            std::cout << "Naive Blocking Blelloch Parallel Time : " << parRuntime << std::endl;
            std::cout << "Naive Blocking Blelloch Speed-Up : " << (serRuntime / parRuntime) << std::endl;
        }
    }

    //Finalise the MPI call
    MPI_Finalize();

    return 0;
}

/*
 * Checks if the provided argument is correct, exits if it isn't
 */
void checkInput(const std::string &arrSize, int threadCount) {
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
 * A function that performs a 1D blocking naive blelloch scan in parallel
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

    //Compute the full scan sum for each thread on its subset of elements
    for (int i = 1; i < localN; i++) {
        localIn[i] += localIn[i - 1];
    }

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

/*
 * A function that performs classic work efficient blelloch on an array in parallel
 */
void mpiBlellochScan(std::vector<int> &in, int N) {
    // Initialise variables for each thread that will be used in the algorithm
    int threadCount, threadID, localN;

    // Get the current thread rank and total thread count
    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCount);

    // Get the size of the local array (size of the array each thread will have)
    localN = N / threadCount;

    // Initialise 2 arrays that will be used to store the full scan sums for each thread
    std::vector<int> localIn(localN, 0);
    std::vector<int> parallelIN(N, 0);

    // Scatter / distribute the main array to each thread into the localIn array (of size localN)
    MPI_Scatter(in.data(), localN, MPI_INT, localIn.data(), localN, MPI_INT, 0, MPI_COMM_WORLD);

    //Set the initial count for communication to 1 (since cant divide by 0)
    int count = 1;

    //Variables needed throughout the up sweep and down sweep
    //d, k are loop variables, inc is an increment variable for the loop
    //ind1 & ind2 are the variables used to index the array
    int d, k, ind1, ind2, i, t, maxNum, shift1, shift2;

    //Up-sweep part of the algorithm
    for (d = 0; d < (int) log2(N); d++) {
        shift1 = 1 << d;

        //Make sure index 1 is within the sub array
        if (shift1 < localN) {
            shift2 = 1 << (d + 1); //Store the incremebt in a variable
            //Use a tree traversal for the local array
            for (k = 0; k < localN - 1; k += shift2) {
                ind2 = k + shift2 - 1;
                ind1 = k + shift1 - 1;
                localIn[ind2] += localIn[ind1];
            }

        } else if (shift1 == localN) { //If the index is the end element in the last array
            //For down sweep, modulus is used to communicate the end values to the next neighour but only for select nodes
            if ((threadID & 1) == 0) {
                MPI_Send(&localIn[localN - 1], 1, MPI_INT, threadID + 1, 0, MPI_COMM_WORLD);
            } else {
                MPI_Recv(&ind2, 1, MPI_INT, threadID - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                localIn[localN - 1] += ind2; //Increment the value by the recieved value
            }
        } else if (shift1 > localN) { //if the index is beyond the current array size
            //Increase the counts needed (this ensures not every thread will
            //communicate but rather only those in the tree traversal)
            count++;
            int tmp = 1 << (count - 1); //Temporary increment variable

            //Modulus is used again to determine which thread to send to (again to conform to that tree structure)
            if ((threadID & (tmp - 1)) == tmp - 1) {
                if (((threadID / tmp) & 1) == 0) {
                    MPI_Send(&localIn[localN - 1], 1, MPI_INT, threadID + tmp, 0, MPI_COMM_WORLD);
                } else {
                    MPI_Recv(&ind2, 1, MPI_INT, threadID - tmp, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //Add on the value to the array index
                    localIn[localN - 1] += ind2;
                }
            }
        }
    }

    //Store the last value from the last local sub array stored in the end most (the highest ranked) thread
    if (threadID == threadCount - 1) {
        maxNum = localIn[localN - 1]; //Backup the end element that will be used for the left shift
        localIn[localN - 1] = 0; //Also set the value for 0 in order to do a pre scan from the tree
    }

    //Proceeding with the down sweep part of the algorithm:

    int sCount = log2(threadCount); //Get the height of the local tree (used for the selective communication)
    //Start the local tree exploration
    for (d = log2(N) - 1; d >= 0; --d) {
        shift1 = 1 << d;
        if (shift1 < localN) {
            shift2 = 1 << (d + 1);
            for (i = 0; i < localN - 1; i += shift2) {
                t = localIn[i + shift1 - 1];
                localIn[i + shift1 - 1] = localIn[i + shift2 - 1];
                localIn[i + shift2 - 1] = t + localIn[i + shift2 - 1];
            }
            //Same as above, checking where in relation to the localN the index is
        } else if (shift1 == localN) {
            int tmp = 1 << (sCount - 1); //Temp variable for the increment
            if ((threadID & (tmp - 1)) == tmp - 1) { //Use modulus to determine the selective communications
                int copyNum = 0; //Define a variable to save the value that gets overwritten in the down sweep when the values get swapped
                int destinationThread = 0;
                //Use modulus again to determine the send and receive pairs
                if (((threadID / tmp) & 1) == 0) {
                    destinationThread = threadID + tmp;
                } else {
                    destinationThread = threadID - tmp;
                }
                //Use a single sendrecv in order to avoid a possible deadlock
                MPI_Sendrecv(&localIn[localN - 1], 1, MPI_INT, destinationThread, 0, &copyNum, 1, MPI_INT,
                             destinationThread, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (((threadID / tmp) & 1) == 0) {
                    localIn[localN - 1] = copyNum; //Repalce the left value of the communcation with the right
                } else {
                    localIn[localN - 1] = copyNum + localIn[localN -
                                                            1]; //Add the left value to the right balue and store in the right
                }
            }
            sCount--;
            //If the index lies outside the local array
        } else if (shift1 > localN) {
            int tmp = 1 << (sCount - 1);
            //Again, use a modulus to denote the spcific tree-like communication
            if ((threadID & (tmp - 1)) == tmp - 1) {
                int copyNum = 0;
                int destinationThread = 0;
                if (((threadID / tmp) & 1) == 0) {
                    destinationThread = threadID + tmp;
                } else {
                    destinationThread = threadID - tmp;
                }
                //Use a send and receive pair to send the elements
                MPI_Sendrecv(&localIn[localN - 1], 1, MPI_INT, destinationThread, 0, &copyNum, 1, MPI_INT,
                             destinationThread, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //Swap the elements around and increase the right element by the left element
                if (((threadID / tmp) & 1) == 0) {
                    localIn[localN - 1] = copyNum;
                } else {
                    localIn[localN - 1] = copyNum + localIn[localN - 1];
                }
            }
            sCount--;
        }
    }

    //Gather  the data back to the main array
    MPI_Allgather(localIn.data(), localN, MPI_INT, parallelIN.data(), localN, MPI_INT, MPI_COMM_WORLD);

    //We then do a scan and a left shift to get the pre scan into an inclusive scan
    if (threadID == threadCount - 1) {
        parallelIN.push_back(maxNum);
        parallelIN.erase(parallelIN.begin());
    }

    //Get the final inclusive scan array back to all threads so that it can be returned for validation
    MPI_Bcast(parallelIN.data(), N, MPI_INT, threadCount - 1, MPI_COMM_WORLD);
    if (threadID == 0) {
        in = parallelIN;
    }


}
