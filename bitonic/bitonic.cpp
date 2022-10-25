#include <iostream>
#include <vector>
#include <random>
#include <math.h>
#include <omp.h>

/*
 * Checks if the provided argument is correct, exits if it isn't
 */
void checkInput(std::string arrSize) {
    if (stoi(arrSize) < 1 || stoi(arrSize) > 28) std::cout << "INCORRECT ARRAY SIZE\n", exit(1);
}

/*
 * Generate an array of size N with random elements between [1, 1000)
 */
std::vector<int> generateArray(int N)
{
    std::vector<int> arr(N, 0);
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<> dist(0, 50);
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
 * Serial Bitonic Sort
 */
void serialBitonic(std::vector<int> &arr)
{
    int N = arr.size();
    // k represents number of times of bitonic converison(merge two bitonic sequence into one), so it is log2(N) times
    for (int k = 0; k < log2(N); ++k)
    {
        // j represents number of times of bitonic-sorting, so it is k times, but j is decreasing
        for (int j = k; j >= 0; --j)
        {
            // each line will run number comparison, all line will run N/2 comparison
            for (int i = 0; i < N / 2; ++i)
            {
                // numberDifference is used for number 1 and number 2 separation and different number1 with the increased on i
                int numDiff = pow(2, j);
                int numTimes = (int)(i / numDiff);
                // get different positions for different i
                int pos = i + numDiff * numTimes;
                // get position for number 1 and number 2 with correct number differnce
                int number1 = pos, number2 = pos + numDiff;

                // multiple just determine whether the number1 should be increasing or decreasing within  the rangeNum
                // if multiple(1) is increasing and multiple(-1) is decreasing
                int rangeNum = pow(2, k);
                int multiple = pow(-1, (int)(i / rangeNum));

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
/*
 *Get times for serial Bitonic
 */
int main(int argc, char *argv[])
{
    //Check the input array size is correct
    checkInput(argv[1]);

    // generate size of the array
    int N = (int)pow(2, atoi(argv[1]));
    std::vector<int> arr = generateArray(N);
    
    // printLine(arr);
    double times = omp_get_wtime();
    serialBitonic(arr);
    times = omp_get_wtime() - times;
    // printLine(arr);
    std::cout << "Serial times: " << times << "\n";
    return 0;
}
