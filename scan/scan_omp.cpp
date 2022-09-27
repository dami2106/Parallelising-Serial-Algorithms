/*
 * Parallel implementation of Blelloch Parallel Scan
 */
#include <iostream>
#include <vector>
#include <array>
#include <random>
//#include <omp.h>

std::vector<int> generateArray(int N);
void fullScan(std::vector<int> &in, std::vector<int> &out, int N);
void ompFullScan(std::vector<int> &in, int N);
bool compareScan(std::vector<int> &arrOne, std::vector<int> &arrTwo, int N);

int main() {
    const int N = 8;
    std::vector<int> in;

    in.push_back(2);
    in.push_back(1);
    in.push_back(4);
    in.push_back(0);
    in.push_back(3);
    in.push_back(7);
    in.push_back(6);
    in.push_back(3);

    for(int c : in)
        std::cout << c << " ";
    std::cout << "\n";
    ompFullScan(in, N);
    for(int c : in)
        std::cout << c << " ";
    std::cout << "\n";
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
void ompFullScan(std::vector<int> &in, int N) {
    for (int d = 0 ; d < (int)log2(N)  ; d++) {
        int inc = (int)pow(2, d+1);
        for(int k = 0 ; k < N-1 ; k += inc) {
            int ind1 = k + inc - 1;
            int ind2 = k + (int)pow(2, d) - 1;

            if(ind1 < N && ind2 < N)
                in[ind1] += in[ind2];

        }
    }


    int temp=in[N-1];
    in[N-1]= 0;


    for (int  d = log2(N)-1 ;  d >=0 ; -- d) {
;
        for (int i = 0; i <=N-1 ; i+= pow(2,d+1)) {
            int t = in[i + pow(2,d)-1] ;
            in[i + pow(2,d)-1] =in[i + pow(2,d+1)-1];
            in[i + pow(2,d+1)-1] = t + in[i + pow(2,d+1)-1];

        }
    }

    in.push_back(temp);
    in.erase(in.begin());
}

/*
 * A function that checks that both arrays are equal (both have the same full scan)
 */
bool compareScan(std::vector<int> &arrOne, std::vector<int> &arrTwo, int N) {
    for(int i = 0 ; i < N ; i++) {
        if(arrOne[i] != arrTwo[i])
            return false;
    }
    return true;
}