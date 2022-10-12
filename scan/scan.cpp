/*
 * Serial implementation of Blelloch Scan
 */
#include <iostream>
#include <vector>
#include <random>
#include <omp.h>

std::vector<int> generateArray(int N);
void fullScan(std::vector<int> &in, std::vector<int> &out, int N);

int main(int argc, char *argv[]) {
    //Get the size of the array from the parameters
    int N = (int) pow(2, atoi(argv[1]));
    //Set the variables used for timing
    double startTime, endTime = 0;

    //Generate an array of random elements of size N
    std::vector<int> in = generateArray(N);

    //An array of 0s of size N used to store the result of full scan
    std::vector<int> out(N, 0);

    //Time and call the serial full scan
    startTime = omp_get_wtime();
    fullScan(in, out, N);
    endTime = omp_get_wtime() - startTime;

    //Return the timing since serial has no speed-up or validation
    std::cout << "Serial Time: " << endTime << std::endl;
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
void fullScan(std::vector<int> &in, std::vector<int> &out, int N) {
    out[0] = in[0]; //Set the first elements to be equal
    //Perform serial full scan
    for (int i = 1; i < N; i++)
        out[i] = in[i] + out[i - 1];
}