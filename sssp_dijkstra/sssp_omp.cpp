#include <iostream>
#include <vector>
#include <sstream>
#include <array>
#include <unordered_set>
#include <omp.h>
#include <fstream>
#include <climits>

using namespace std;

#define NUMTHREADS 8

vector<vector<int> > makeGraph(int &vertexCount, int &edgeCount, const string &fileName);

vector<int> serialDijkstra(int vertexCount, int startVertex, vector<vector<int> > adj);

vector<int> parallelDijkstra(int vertexCount, int startVertex, vector<vector<int> > adj);

void printPath(int vert, vector<int> parents);

void printSolution(int startVertex, vector<int> distances, vector<int> parents);

//arg[1] = file name, arg[2] = iterations
int main(int argc, char *argv[]) {
    int vertexCount, edgeCount, startVertex = 0;

    vector<vector<int> > adj = makeGraph(vertexCount, edgeCount, argv[1]);
    double startTime, serRunTime = 0, parRunTime = 0;

    startTime = omp_get_wtime();
    vector<int> serialDist = serialDijkstra(vertexCount, startVertex, adj);
    serRunTime += omp_get_wtime() - startTime;

    startTime = omp_get_wtime();
    vector<int> parallelDist = parallelDijkstra(vertexCount, startVertex, adj);
    parRunTime += omp_get_wtime() - startTime;

    for (auto c: serialDist)
        cout << c << " ";
    cout << endl;
    for (auto c: parallelDist)
        cout << c << " ";
    cout << endl;

    if (serialDist != parallelDist)
        cout << "(Validation Failed!)";
    else
        cout << serRunTime / parRunTime << "  (Validation Passed!)";
//    cout << "Serial Time : " << serRunTime / iterations << endl;
//    cout << "Parallel Time : " << parRunTime / iterations << endl;
//    cout << "Speed-Up : " << (serRunTime / parRunTime) << endl;
}

vector<int> serialDijkstra(int vertexCount, int startVertex, vector<vector<int> > adj) {
    unordered_set<int> vT;
    vector<int> l(vertexCount, INT_MAX);
    //vector<int> parents(vertexCount);
    l[startVertex] = 0;
    //parents[startVertex] = -1;

    while ((int) vT.size() != vertexCount) {
        int u, min = INT_MAX;

        for (int i = 0; i < vertexCount; i++)
            if (vT.find(i) == vT.end() && l[i] < min) {
                min = l[i];
                u = i;
            }

        vT.insert(u);

        for (int v = 0; v < vertexCount; v++) {
            if (adj[v][u] != INT_MAX)
                if (vT.find(v) == vT.end() && l[v] > l[u] + adj[v][u]) {
                    l[v] = l[u] + adj[v][u];
                    //parents[v] = u;
                }
        }
    }
    return l;
}

vector<int> parallelDijkstra(int vertexCount, int startVertex, vector<vector<int> > adj) {
    //unordered_set<int> vT; //Keeps track of vertices explored
    vector<bool> vTn(vertexCount, false);
    vector<int> l(vertexCount, INT_MAX); //Need to subset per thread - keeps track of the min distance from startVertex

    int threadID, threadCount, localMin = INT_MAX, localU = -1, currentVert, threadBoundLeft, threadBoundRight;
    int u = -1, min = INT_MAX;

    l[startVertex] = 0;

#pragma omp parallel num_threads(NUMTHREADS) firstprivate(localMin, localU) private(threadID, threadCount, currentVert, threadBoundLeft, threadBoundRight) shared(vTn, min, u, adj, l, startVertex)
    {
        threadID = omp_get_thread_num();
        threadCount = omp_get_num_threads();
        threadBoundLeft = threadID * (vertexCount / threadCount);
        threadBoundRight = ((threadID + 1) * (vertexCount / threadCount)) - 1;

        for (currentVert = 0; currentVert < vertexCount; currentVert++) {
#pragma omp single
            {
                u = -1;
                min = INT_MAX;
            }

            localU = -1;
            localMin = INT_MAX;

            for (int i = threadBoundLeft; i <= threadBoundRight; i++)
                if (!vTn[i] && l[i] < localMin) localMin = l[i], localU = i;

#pragma omp critical
            {
                if (localMin < min) {
                    min = localMin;
                    u = localU;
                }
            }
#pragma omp barrier

#pragma omp single
            if (u != -1) vTn[u] = true;
#pragma omp barrier


            if (u != -1) {
                for (int i = threadBoundLeft; i <= threadBoundRight; i++) {
                    if (adj[i][u] != INT_MAX)
                        if (!vTn[i] && l[i] > l[u] + adj[i][u]) l[i] = l[u] + adj[i][u];
                }
            }

#pragma omp barrier
        }
    }
    return l;
}


/*
 * Prints the path from start to vertex
 */
void printPath(int vert, vector<int> parents) {
    if (vert == -1)
        return;
    printPath(parents[vert], parents);
    cout << vert << " ";
}

/*
 * Prints the entirety of the solution
 */
void printSolution(int startVertex, vector<int> distances, vector<int> parents) {
    int nVertices = distances.size();

    cout << "The distance from " << startVertex << " to each vertex is:\n";
    cout << "  v    dist " << startVertex << "->v\n----   ---------\n";
    for (int vertexIndex = 0; vertexIndex < nVertices; vertexIndex++)
        if (vertexIndex != startVertex) cout << vertexIndex << ":\t " << distances[vertexIndex] << "\n";

    cout << "\n\nThe shortest path from " << startVertex << " to each vertex is:\n";
    cout << "  v    dist " << startVertex << "->v\n----   ---------\n";
    for (int vertexIndex = 0; vertexIndex < nVertices; vertexIndex++) {
        if (vertexIndex != startVertex) {
            cout << vertexIndex << ":\t ";
            printPath(vertexIndex, parents);
            cout << "\n";
        }
    }
}

/*
 * Creates and returns a square symmetric adjacency matrix from file in
 */
vector<vector<int> > makeGraph(int &vertexCount, int &edgeCount, const string &fileName) {
    string initLine, word;
    vector<int> graphInfo;

    ifstream fileReader("graphs/" + fileName + ".txt");

    getline(fileReader, initLine);
    stringstream ss(initLine);

    while (ss >> word)
        graphInfo.push_back(stoi(word));

    vertexCount = graphInfo[0];
    edgeCount = graphInfo[1];

    vector<vector<int> > adj(vertexCount, vector<int>(vertexCount, INT_MAX));
    for (int i = 0; i < edgeCount; ++i) {
        string currNode, nodeItem;
        array<int, 3> nodeInfo;
        int k = 0;

        getline(fileReader, currNode);
        stringstream sx(currNode);

        while (sx >> nodeItem) {
            nodeInfo[k] = stoi(nodeItem);
            ++k;
        }
        adj[nodeInfo[0]][nodeInfo[1]] = nodeInfo[2];
        adj[nodeInfo[1]][nodeInfo[0]] = nodeInfo[2];
    }

    fileReader.close();
    return adj;
}
