/*
 * Parallel implementation of Blelloch Parallel Scan
 */
#include <iostream>
#include <vector>
#include <random>
#include <omp.h>

std::vector<int> generateArray(int N);

void fullScan(std::vector<int> &in, std::vector<int> &out, int N);

void ompFullScan(std::vector<int> &in, std::vector<int> &out, int N);

int main(int argc, char *argv[]) {

    int N = atoi(argv[1]);
    int iterations = atoi(argv[2]);

    double startTime, runTime = 0;
//    for(int iter = 0 ; iter < iterations ; iter++) {
//        std::vector<int> in = generateArray(N);
//        std::vector<int> out(N, 0);
//
//        startTime = omp_get_wtime();
//        fullScan(in, out, N);
//        runTime += omp_get_wtime() - startTime;
//    }
    std::vector<int> in = generateArray(N);
    std::vector<int> out(N, 0);
    for (int c: in)
        std::cout << c << " ";
    std::cout << "\n";
    ompFullScan(in, out, N);
    for (int c: in)
        std::cout << c << " ";


    std::cout << "Serial time: " << (runTime / iterations);
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

/*
 * A function that performs a parallel full scan on the given array
 */
void ompFullScan(std::vector<int> &in, std::vector<int> &out, int N) {
    int twoip1 = 0, twoi = 0;

    for (int i = 0; i <= log2(N); i++) {
        twoip1 = 1 << (i + 1);
        twoi = 1 << i;

        for (int j = 0; j < N; j += twoip1) {
            in[j + twoip1 - 1] = in[j + twoi - 1] + in[j + twoip1 - 1];
        }
    }



//    for (int i = (int)log2(N); i >= 0; i--) {
//        twoip1 = 1 << (i + 1);
//        twoi = 1 << i;
//
//        for (int j = 0; j < N; j += twoip1) {
//            long t = out[j + twoi - 1];
//            out[j + twoi - 1] = out[j + twoip1 - 1];
//            out[j + twoip1 - 1] = t + out[j + twoip1 - 1];
//        }
//    }
}