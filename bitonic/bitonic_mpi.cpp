#include <iostream>
#include <vector>
#include <random>
#include <math.h>
#include <mpi.h>


/*
 * Checks if the provided argument is correct, exits if it isn't
 */
void checkInput(std::string arrSize, int threadCount) {
    if (pow(2, stoi(arrSize)) < threadCount) std::cout << "INCORRECT ARRAY SIZE\n", exit(1);
    if (stoi(arrSize) < 3 || stoi(arrSize) > 28) std::cout << "INCORRECT ARRAY SIZE\n", exit(1);
}

 /*
  * Generate an array of size N with random elements between [1, 1000)
  */
 std::vector<int> generateArray(int N)
 {
     std::vector<int> arr(N, 0);
     std::random_device rd;
     std::mt19937 generator(rd());
     std::uniform_int_distribution<> dist(0, 1000);
     for (int i = 0; i < N; i++)
         arr[i] = dist(generator);
     return arr;
 }


/*
 * swap the number at position 1 and position 2 around
 */
void swapNum(std::vector<int> &arr, int pos1, int pos2)
{
    int temp = arr[pos1];
    arr[pos1] = arr[pos2];
    arr[pos2] = temp;
}
/*
 *print out all vector numbers
 */
void printLine(std::vector<int> &arr)
{
    for (int c : arr)
        std::cout << c << " ";
    std::cout << "\n";
}
/*
 * Serial Bitonic
 */
void serialBitonic(std::vector<int> &arr)
{
    int N = arr.size();
    // k represents number of times of bitonic merge, so it is log2(N) times
    for (int k = 0; k < log2(N); ++k)
    {
        // j represents number of times of bitonic sorts, so it is k times, but j is decreasing
        for (int j = k; j >= 0; --j)
        {
            // each line will run number comparison, all line will run N/2 comparison
            for (int i = 0; i < N / 2; ++i)
            {
                // numberDifference is used for number 1 and number 2 separation and different number1 with the increased on i
                int numDiff = pow(2, j);
                int numTimes = (int)(i / numDiff);

                // get different position for different i
                int pos = i + numDiff * numTimes;
                // get position for number 1 and number 2 with correct number differnce
                int number1 = pos, number2 = pos + numDiff;

                // multiple just determine whether the comparison between two number should be increasing or decreasing
                // if multiple(1),it is increasing and multiple(-1), it is decreasing
                int multiple = pow(-1, (int)(i / pow(2, k)));

                // two scenarios:
                // with multiple(1), number 1 is greater than number 2 , we swap number
                // with multiple(-1), number 1 is greater than number 2(actually without multiple, number 1 is smaller than number2) , we swap number
                if (arr[number1] * (multiple) > arr[number2] * multiple)
                {
                    swapNum(arr, number2, number1);
                }
            }
        }
    }
}
/**
 * parallel mpi bitonic implementation
 */
void parallelBitonic(std::vector<int> &arr)
{
    int threadCount, threadID;
    double parallelTime = 0, serialTime = 0;
    int N = arr.size();
    // Initialize MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCount);
    MPI_Status status;

    // Save a copy of the vector into serialArr, for serial bitonic
    std::vector<int> seiralArr = arr;
    
    // divideNUmbers represents the numbers of element get in each thread
    int divideNumbers = N / threadCount;
    // calculated parallel bitonic timing and start timer for MPI
    if (threadID == 0)
    {
        serialTime = MPI_Wtime();
        serialBitonic(seiralArr);
        serialTime = MPI_Wtime() - serialTime;
        parallelTime = MPI_Wtime();
    }

    // create a temporary vector with numbers required in each threads
    std::vector<int> mpiArr(divideNumbers, 0);
    // Distribute the arr data equally to mpiArr, so each thread will get divideNumbers of elements
    MPI_Scatter(arr.data(), divideNumbers, MPI_INT, mpiArr.data(), divideNumbers, MPI_INT, 0, MPI_COMM_WORLD);

    // k represents number of times of bitonic merge, so it is log2(N) times
    for (int k = 0; k < log2(N) ; ++k)
    {
        // j represents number of times of bitonic sorts, so it is k times, but j is decreasing
        for (int j = k; j >= 0; --j)
        {
            // if elements in the thread requires another elements in other thread(mpi sendRec)
            if ((2 * pow(2, j)) > divideNumbers)
            {   
                //initaliszed a copy array to receive the data from another threads
                std::vector<int> copyArr(divideNumbers, 0);
                //division is thread number difference and determine which division is current thread on
                int division = pow(2,j)/divideNumbers;
                int destinationThread = 0;
                //decision is determine whether current thread has to communicate forward or backward
                int decision = (int)threadID / division;
                if (decision % 2 == 0)
                {
                    destinationThread = threadID + division;
                }
                else
                {
                    destinationThread = threadID - division;
                }
                
                // send mpiArr data to destinationThread and receive mpiArr data to copyArr from destinationThread
                MPI_Sendrecv(mpiArr.data(), divideNumbers, MPI_INT, destinationThread, 0, copyArr.data(), divideNumbers, MPI_INT, destinationThread, 0, MPI_COMM_WORLD, &status);

                //get multiple of both thread whether they are increasing or decreasing
                int multiple = pow(-1, (int)((threadID*divideNumbers) / (pow(2,k)*2)));
                //if cuurent thread number is less than destination thread number, we mutiply multiple by -1
                if (threadID < destinationThread)
                {
                    multiple = multiple * -1;
                }
                //comparing mpiArr data to copyArr, and decide whether to swap
                for (int i = 0; i < divideNumbers; i++)
                {
                    if (mpiArr[i] * (multiple) < copyArr[i] * multiple)
                    {
                        int temp = copyArr[i];
                        mpiArr[i] = temp;
                    }
                }
            }
            else
            {
                // N/2 comparsion will run on each threads 
                for (int i = 0; i < divideNumbers / 2; i++)
                {
                    // numberDifference is used for number 1 and number 2 separation and different number1 with the increased on i
                    int numDiff = pow(2, j);
                    int numTimes = (int)(i / numDiff);

                    // get different position for different i
                    int pos = i + numDiff * numTimes;
                    // get position for number 1 and number 2 with correct number differnce
                    int number1 = pos, number2 = pos + numDiff;

                    // multiple just determine whether the comparison between two number should be
                    //  increasing or decreasing by determine the oringinal position in the array
                    int rangeNum = pow(2, k)*2;
                    int multiple = pow(-1, (int)((pos + (divideNumbers) * threadID) / rangeNum));

                    // two scenarios:
                    // with multiple(1), number 1 is greater than number 2 , we swap number
                    // with multiple(-1), number 1 is greater than number 2(actually without multiple, number 1 is smaller than number2) , we swap number
                    if (mpiArr[number1] * (multiple) > mpiArr[number2] * multiple)
                    {
                        swapNum(mpiArr, number2, number1);
                    }
                }
            }
        }
    }
    // group all the data back into arr
    MPI_Allgather(mpiArr.data(), divideNumbers, MPI_INT, arr.data(), divideNumbers, MPI_INT, MPI_COMM_WORLD);
    // verify is the parallel arr is the same as serialArr that the answer is correct


    if (threadID == 0) {
        parallelTime = MPI_Wtime() - parallelTime;

        //Validate the parallel data against the serial sum array and serial parallel array
        if (seiralArr != arr) {
            std::cout << "(Validation Failed!)\n";
        } else {
            std::cout << "(Validation Passed!)\n";
            std::cout << "Serial Time : " << serialTime << std::endl;
            std::cout << "Parallel Time : " << parallelTime << std::endl;
            std::cout << "Speed-Up : " << (serialTime / parallelTime) << std::endl;
            std::cout << "Efficiency : " << (serialTime / parallelTime)/threadCount << std::endl;
        }
    }

    // Close the MPI
    MPI_Finalize();
}
int main(int argc, char *argv[])
{

    //Check the input array size is correct
    //checkInput(argv[1], threadCount);

    // generate arr
    int N = (int)pow(2, atoi(argv[1]));

    //Serial + time
    std::vector<int> arr = generateArray(N);
    
    // run parallel bitonic
    parallelBitonic(arr);

    return 0;
}