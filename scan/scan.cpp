/*
 * Serial implementation of Blelloch Parallel Scan
 */
#include <iostream>
#include <vector>
#include <random>
#include <omp.h>

std::vector<int> generateArray(int N);
void fullScan(std::vector<int> &in, std::vector<int> &out, int N);

int main(int argc, char *argv[]) {
    int N = (int) pow(2, atoi(argv[1]));
    double startTime, endTime = 0;

    std::vector<int> in = generateArray(N);

    std::vector<int> out(N, 0);

    startTime = omp_get_wtime();
    fullScan(in, out, N);
    endTime = omp_get_wtime() - startTime;

    std::cout << (endTime);
    return 0;
}

/*
 * Generate an array of size N with random elements between [1, 50)
 */
std::vector<int> generateArray(int N) {
    std::vector<int> arr(N, 0);

    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<> dist(1, 50);

    for (int i = 0; i < N; i++)
        arr[i] = dist(generator);

    return arr;
}

/*
 * A function that performs a serial full scan on the given array
 */
void fullScan(std::vector<int> &in, std::vector<int> &out, int N) {
    out[0] = in[0];
    for (int i = 1; i < N; i++)
        out[i] = in[i] + out[i - 1];
}