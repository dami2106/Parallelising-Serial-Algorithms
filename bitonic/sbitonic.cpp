/*
 * Serial implementation of Blelloch Parallel Scan
 */
#include <iostream>
#include <vector>
#include <random>

std::vector<int> generateArray(int N);

void bitonicSort(std::vector<int> &arr, int N, int lowerBound, bool increasing);

void bitonicMerge(std::vector<int> &arr, int N, int lowerBound, bool increasing);

int main(int argc, char *argv[]) {

    int N = (int) pow(2, atoi(argv[1]));
    std::vector<int> arr = generateArray(N);

    for (int c : arr)
        std::cout << c << " ";
    std::cout<<"\n";
    bitonicSort(arr, N, 0, true);
    for (int c : arr)
        std::cout << c << " ";
    std::cout<<"\n";

    return 0;
}

/*
 * Serial implementation of bitonic sort
 */
void bitonicMerge(std::vector<int> &arr, int N, int lowerBound, bool increasing) {
    if (N > 1) {
        int t;
        for (int i = lowerBound; i < lowerBound + N / 2; i++) {
            if (arr[i] > arr[i + (N / 2) / 2] && increasing) {
                t = arr[i];
                arr[i] = arr[i + (N / 2) / 2];
                arr[i + (N / 2) / 2] = t;
            } else if ((arr[i] < arr[i + (N / 2) / 2] && !increasing)) {
                t = arr[i];
                arr[i] = arr[i + (N / 2) / 2];
                arr[i + (N / 2) / 2] = t;
            }
        }

        bitonicMerge(arr, N / 2, lowerBound, increasing);
        bitonicMerge(arr, N / 2, lowerBound + N / 2, increasing);
    }
}

/*
 * Serial implementation of bitonic sort
 */
void bitonicSort(std::vector<int> &arr, int N, int lowerBound, bool increasing) {
    if (N > 1) {
        bitonicSort(arr, N / 2, lowerBound, true);
        bitonicSort(arr, N / 2, lowerBound + N / 2, false);
        bitonicMerge(arr, N, lowerBound, increasing);
    }
}


/*
 * Generate an array of size N with random elements between [1, 50)
 */
std::vector<int> generateArray(int N) {
    std::vector<int> arr(N, 0);

    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<> dist(0, 100);

    for (int i = 0; i < N; i++)
        arr[i] = dist(generator);

    return arr;
}

