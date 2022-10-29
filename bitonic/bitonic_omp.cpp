#include <iostream>
#include <vector>
#include <random>
#include <math.h>
#include <omp.h>

#define NUMTHREADS 8

/*
 * Generate an array of size N with random elements between [1, 1000)
 */
std::vector<int> generateArray(int N) {
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
void swapNum(std::vector<int> &arr, int pos1, int pos2) {
    int temp = arr[pos1];
    arr[pos1] = arr[pos2];
    arr[pos2] = temp;
}

/*
 *print out all vector numbers
 */
void printLine(std::vector<int> &arr) {
    for (int c: arr)
        std::cout << c << " ";
    std::cout << "\n";
}

/*
 * Serial Bitonic
 */
void serialBitonic(std::vector<int> &arr) {
    int N = arr.size();
    // k represents number of times of bitonic converison(merge two bitonic sequence into one), so it is log2(N) times
    for (int k = 0; k < log2(N); ++k) {
        // j represents number of times of bitonic-sorting, so it is k times, but j is decreasing
        for (int j = k; j >= 0; --j) {
            // each line will run number comparison, all line will run N/2 comparison
            for (int i = 0; i < N / 2; ++i) {
                // numberDifference is used for number 1 and number 2 separation and different number1 with the increased on i
                int numDiff = pow(2, j);
                int numTimes = (int) (i / numDiff);

                // get different position for different i
                int pos = i + numDiff * numTimes;
                // get position for number 1 and number 2 with correct number differnce
                int number1 = pos, number2 = pos + numDiff;

                // multiple just determine whether the comparison between two number should be increasing or decreasing 
                // if multiple(1),it is increasing and multiple(-1), it is decreasing
                int multiple = pow(-1, (int) (i / pow(2, k)));

                // two scenarios:
                // with multiple(1), number 1 is greater than number 2 , we swap number
                // with multiple(-1), number 1 is greater than number 2(actually without multiple, number 1 is smaller than number2) , we swap number
                if (arr[number1] * (multiple) > arr[number2] * multiple) {
                    swapNum(arr, number2, number1);
                }
            }
        }
    }
}

/*
 * Parallel Bitonic
 */
void parallelBitonic(std::vector<int> &arr) {
    int N = arr.size();
    // k represents number of times of bitonic converison(merge two bitonic sequence into one), so it is log2(N) times
    for (int k = 0; k < log2(N); ++k) {
        // j represents number of times of bitonic-sorting, so it is k times, but j is decreasing
        for (int j = k; j >= 0; --j) {
            // using parallel for to parallize the follwing for loop, so the N/2 comparsion will run in parallel
#pragma omp parallel for num_threads(NUMTHREADS) shared(N, arr, k, j) default(none)
            for (int i = 0; i < N / 2; ++i) {
                // numberDifference is used for number 1 and number 2 separation and different number1 with the increased on i
                int numDiff = pow(2, j);
                int numTimes = (int) (i / numDiff);

                // get different position for different i
                int pos = i + numDiff * numTimes;
                // get position for number 1 and number 2 with correct number differnce
                int number1 = pos, number2 = pos + numDiff;

                // multiple just determine whether the comparison between two number should be increasing or decreasing 
                // if multiple(1),it is increasing and multiple(-1), it is decreasing
                int multiple = pow(-1, (int) (i / pow(2, k)));

                // two scenarios:
                // with multiple(1), number 1 is greater than number 2 , we swap number
                // with multiple(-1), number 1 is greater than number 2(actually without multiple, number 1 is smaller than number2) , we swap number
                if (arr[number1] * (multiple) > arr[number2] * multiple) {
                    swapNum(arr, number2, number1);
                }
            }
        }
    }
}

/*
 * Checks if the provided argument is correct, exits if it isn't
 */
void checkInput(std::string arrSize) {
    if (pow(2, stoi(arrSize)) < NUMTHREADS) std::cout << "INCORRECT ARRAY SIZE\n", exit(1);
    if (stoi(arrSize) < 3 || stoi(arrSize) > 28) std::cout << "INCORRECT ARRAY SIZE\n", exit(1);
}

int main(int argc, char *argv[]) {

    //Check the input array size is correct
    checkInput(argv[1]);

    // generate size of the array
    int N = (int) pow(2, atoi(argv[1]));
    std::vector<int> arr = generateArray(N);

    //calculate time for serial implementation
    std::vector<int> seriesArr = arr;
    double sTimes = omp_get_wtime();
    serialBitonic(seriesArr);
    sTimes = omp_get_wtime() - sTimes;


    //calculate time for parallel implementation
    std::vector<int> parallelArr = arr;
    double pTimes = omp_get_wtime();
    parallelBitonic(parallelArr);
    pTimes = omp_get_wtime() - pTimes;

    // verify to see is the answer correct
    if (parallelArr != seriesArr) {
        std::cout << "(Validation Failed!)\n";
    } else {
        std::cout << "(Validation Passed!)\n";
        std::cout << "Serial Time : " << sTimes << std::endl;
        std::cout << "Parallel Time : " << pTimes << std::endl;
        std::cout << "Speed-Up : " << (sTimes / pTimes) << std::endl;
        std::cout << "Efficiency : " << (sTimes / pTimes) / NUMTHREADS << std::endl;
    }

    return 0;
}
