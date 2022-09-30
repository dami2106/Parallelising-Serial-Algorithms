/*
 * Serial implementation of Blelloch Parallel Scan
 */
#include <iostream>
#include <vector>
#include <random>

std::vector<int> generateArray(int N);

int main (int argc, char *argv[]) {

    int N = (int)pow(2, atoi(argv[1]));
    std::vector<int> arr = generateArray(N);

    return 0;
}
void bitonicSort(std::vector<int> &arr, int N, int bottom, int dir){
    if(N > 1) {
        int half = N/2;
        bitonicSort(arr, half, bottom, 1);
        bitonicSort(arr, half, bottom + half, 0);
    }
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

