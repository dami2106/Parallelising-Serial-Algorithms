/*
 * Parallel implementation of Blelloch Parallel Scan
 */
#include <iostream>
#include <vector>
#include <random>
#include <omp.h>

#define NUMTHREADS 8

std::vector<int> generateArray(int N);

void fullScan(std::vector<int> &in, std::vector<int> &out, int N);

void ompFullScan(std::vector<int> &in, int N);

bool compareScan(std::vector<int> &arrOne, std::vector<int> &arrTwo, int N);

int main(int argc, char *argv[]) {
    int N = (int) pow(2, atoi(argv[1]));
    int iter = atoi(argv[2]);
    double startTime, sRunTime = 0, pRunTime = 0;

    std::vector<int> in = generateArray(N);
    for(int c : in)
        std::cout << c << " ";
    std::cout << std::endl;
    ompFullScan(in, N);
    for(int c : in)
        std::cout << c << " ";
    std::cout << std::endl;
//    for (int z = 0; z < iter; z++) {
//        std::vector<int> in = generateArray(N);
//        std::vector<int> out(N, 0);
//
////        startTime = omp_get_wtime();
////        fullScan(in, out, N);
////        sRunTime += omp_get_wtime() - startTime;
//
//        startTime = omp_get_wtime();
//        ompFullScan(in, N);
//        pRunTime += omp_get_wtime() - startTime;
//
//        if (!compareScan(in, out, N))
//            std::cout << "WRONG\n";
//    }


//    std::cout << "Parallel FS gets: " << pRunTime / iter << "\nSerial FS gets: " << sRunTime / iter
//              << "\nWith a speed-up of: " << sRunTime / pRunTime << std::endl;
//    std::cout << iter << " iterations used, for a list of size: 2^" << argv[1] << std::endl;
//    std::cout << "Running on " << NUMTHREADS << " threads\n";

    return 0;
}

/*
 * Generate an array of size N with random elements between [1, 50)
 */
std::vector<int> generateArray(int N) {
    std::vector<int> arr(N, 0);

    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<> dist(0, 10);

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
void ompFullScanOld(std::vector<int> &in, int N) {
    for (int d = 0; d < log2(N); d++) {
        int inc = (int) pow(2, d + 1);

#pragma omp parallel for
        for (int k = 0; k < N; k += inc) {
            in[k + inc - 1] += in[k + (int) pow(2, d) - 1];
        }
    }

    int end = in[N - 1];
    in[N - 1] = 0;


    for (int d = log2(N) - 1; d >= 0; --d) {
        int inc = (int) pow(2, d + 1);
#pragma omp parallel for
        for (int i = 0; i < N; i += inc) {
            int t = in[i + pow(2, d) - 1];
            in[i + pow(2, d) - 1] = in[i + pow(2, d + 1) - 1];
            in[i + pow(2, d + 1) - 1] += t;
        }
    }

    in.push_back(end);
    in.erase(in.begin());
}

/*
 * A function that performs a parallel full scan on the given array
 */
void ompFullScan(std::vector<int> &in, int N) {
    int threadID, threadCount, threadBoundLeft, threadBoundRight, i;
    std::vector<int> globalSum(N, 0);
    std::vector<int> incrementValues(N, 0);

#pragma omp parallel private(i, threadCount, threadID, threadBoundLeft, threadBoundRight) shared(in, N, globalSum, incrementValues)
    {
        threadID = omp_get_thread_num();
        threadCount = omp_get_num_threads();
        threadBoundLeft = threadID * (N / threadCount);
        threadBoundRight = ((threadID + 1) * (N / threadCount)) - 1;

        for(i = threadBoundLeft ; i <= threadBoundRight ; i++)
            in[i] += in[i-1];
        globalSum[threadID] = in[i-1];

#pragma omp barrier

        for(i = 1 ; i < NUMTHREADS ; i++){
            if(threadID >= i){
                incrementValues[threadID] = globalSum[threadID] + globalSum[threadID - 1];
            }
        }
#pragma omp barrier
        for(i = threadBoundLeft ; i <= threadBoundRight ; i++){
            in[i] += globalSum[threadID] - in[threadBoundRight - 1];
        }

    }
}

/*
 * A function that checks that both arrays are equal (both have the same full scan)
 */
bool compareScan(std::vector<int> &arrOne, std::vector<int> &arrTwo, int N) {
    for (int i = 0; i < N; i++) {
        if (arrOne[i] != arrTwo[i])
            return false;
    }
    return true;
}