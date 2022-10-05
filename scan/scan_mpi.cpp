/*
 * Parallel distributed implementation of Blelloch Parallel Scan
 */
#include <iostream>
#include <vector>
#include <random>
#include <mpi.h>

std::vector<int> generateArray(int N);

void fullScan(std::vector<int> &in, std::vector<int> &out, int N);

void mpiFullScan(std::vector<int> &in, int N);

bool compareScan(std::vector<int> &arrOne, std::vector<int> &arrTwo, int N);

int main(int argc, char *argv[]) {
    int N = (int) pow(2, atoi(argv[1]));
    int id;
    double startTime, serRuntime = 0, parRuntime = 0;
    std::vector<int> in;
    std::vector<int> ser(N, 0);

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if (id == 0) {
        in = generateArray(N);
        startTime = MPI_Wtime();
        fullScan(in, ser, N);
        serRuntime = MPI_Wtime() - startTime;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(id == 0)
        startTime = MPI_Wtime();

    mpiFullScan(in, N);

    if(id == 0)
        parRuntime = MPI_Wtime() - startTime;

    MPI_Finalize();

    if (id == 0) {
        if (!compareScan(in, ser, N))
            std::cout << "(Validation Failed!)";
        else {
            std::cout << serRuntime / parRuntime << "  (Validation Passed!)";
//            std::cout << "Parallel FS gets: " << parRuntime / iter << "\nSerial FS gets: " << serRuntime / iter
//                      << "\nWith a speed-up of: " << serRuntime / parRuntime << std::endl;
//            std::cout << iter << " iterations used, for a list of size: 2^" << argv[1] << std::endl;
        }
    }


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
 * A function that performs a distributed parallel full scan on the given array
 */
void mpiFullScan(std::vector<int> &in, int N) {
    int threadCount, threadID, localN, localSum = 0, localIncrement = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCount);
    localN = N / threadCount;

    std::vector<int> globalSum(threadCount, 0);
    std::vector<int> localIn(localN, 0);

    MPI_Scatter(in.data(), localN, MPI_INT, localIn.data(), localN, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 1; i < localN; i++) {
        localIn[i] += localIn[i - 1];
    }
    localSum = localIn[localN - 1];

    //MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(&localSum, 1, MPI_INT, globalSum.data(), 1, MPI_INT, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);


    //Need to sum up the globalSums
    for (int i = 0; i < threadID; i++) {
        localIncrement += globalSum[i];
    }

    if (localIncrement > 0) {
        for (int i = 0; i < localN; i++)
            localIn[i] += localIncrement;
    }

    //MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(localIn.data(), localN, MPI_INT, in.data(), localN, MPI_INT, 0, MPI_COMM_WORLD);
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