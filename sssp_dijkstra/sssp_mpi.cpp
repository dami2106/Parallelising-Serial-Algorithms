#include <iostream>
#include <vector>
#include <sstream>
#include <array>
#include <unordered_set>
#include <fstream>
#include <climits>
#include <mpi.h>
#include <omp.h>

using namespace std;

#define START 0

vector<vector<int>> makeGraph(int &vertexCount, int &edgeCount, const string &fileName);
void serialDijkstra(int vertexCount, vector<vector<int> > adj, vector<int> &l, vector<int> &parents);
void parallelDijkstra(int vertexCount, vector<int> localAdj, vector<int> &globalDist, vector<int> &globalParent);
void printPath(int vert, vector<int> parents);
void printSolution(int startVertex, vector<int> distances, vector<int> parents);


int main(int argc, char *argv[]) {
    double startTime = 0, parRunTime = 0, serRunTime = 0;
    int vertexCount, edgeCount, threadID, threadCount, flatSize, flatDist;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCount);

    vector<int> serDij;
    vector<int> flatAdj;
    vector<vector<int>> adj;
    if (threadID == 0) {
        adj = makeGraph(vertexCount, edgeCount, argv[1]);

        flatSize = vertexCount * vertexCount;
        flatDist = flatSize / threadCount;

        vector<int> temp(flatSize, INT_MAX);
        flatAdj = temp;

        for (int i = 0; i < vertexCount; i++) {
            for (int j = 0; j < vertexCount; j++) {
                flatAdj[i * vertexCount + j] = adj[i][j];
            }
        }
    }

    MPI_Bcast(&flatSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&flatDist, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&vertexCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<int> localAdj(flatDist, INT_MAX);
    vector<int> parallelDist(vertexCount, INT_MAX);
    vector<int> parallelParent(vertexCount, 0);

    vector<int> serialDist(vertexCount, INT_MAX);
    vector<int> serialParent(vertexCount, 0);

    MPI_Scatter(flatAdj.data(), flatDist, MPI_INT, localAdj.data(), flatDist, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (threadID == 0) startTime = omp_get_wtime();
    if (threadID == 0) serialDijkstra(vertexCount, adj, serialDist, serialParent);
    if (threadID == 0) serRunTime = omp_get_wtime() - startTime;

    MPI_Barrier(MPI_COMM_WORLD);

    if (threadID == 0) startTime = omp_get_wtime();
    parallelDijkstra(vertexCount, localAdj, parallelDist, parallelParent);
    if (threadID == 0) parRunTime = omp_get_wtime() - startTime;

    MPI_Barrier(MPI_COMM_WORLD);


    if (threadID == 0) {
        if((serialDist != parallelDist) || (serialParent != parallelParent))
            cout << "(Validation Failed!)";
        else
            cout << serRunTime/parRunTime << "  (Validation Passed!)";
    }

    MPI_Finalize();
}

void serialDijkstra(int vertexCount, vector<vector<int> > adj, vector<int> &l, vector<int> &parents) {
    unordered_set<int> vT;
    l[START] = 0;

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
                    parents[v] = u;
                }
        }
    }
}

void parallelDijkstra(int vertexCount, vector<int> localAdj, vector<int> &globalDist, vector<int> &globalParent) {

    int threadCount, threadID;
    int localVals[2] = {INT_MAX, -1};
    int globalVals[2] = {INT_MAX, -1};

    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCount);

    int localCount = vertexCount / threadCount;
    int lowerBound = threadID * localCount;
    int upperBound = ((threadID + 1) * localCount) - 1;

    vector<bool> visited(localCount, false);
    vector<int> dist(localCount, INT_MAX);
    vector<int> parent(localCount, -1);

    if (lowerBound <= START && START <= upperBound) {
        dist[START] = 0;
        parent[START] = 0;
    }

    for (int i = 0; i < vertexCount; i++) {
        localVals[0] = INT_MAX;
        localVals[1] = -1;

        for (int j = 0; j < localCount; j++) {
            if (!visited[j] && (dist[j] < localVals[0])) {
                localVals[1] = j + lowerBound;
                localVals[0] = dist[j];
            }
        }
        MPI_Allreduce(localVals, globalVals, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);

        if (globalVals[1] == localVals[1]) {
            visited[(globalVals[1] - lowerBound)] = true;
        }


        for (int j = 0; j < localCount; j++) {
            int t = localAdj[j * vertexCount + globalVals[1]];
            if (t != INT_MAX) {
                if (!visited[j] && ((globalVals[0] + t) < dist[j])) {
                    dist[j] = globalVals[0] + t;
                    parent[j] = globalVals[1];
                }
            }

        }

    }

    MPI_Gather(dist.data(), localCount, MPI_INT, globalDist.data(), localCount, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(parent.data(), localCount, MPI_INT, globalParent.data(), localCount, MPI_INT, 0, MPI_COMM_WORLD);
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
vector<vector<int>> makeGraph(int &vertexCount, int &edgeCount, const string &fileName) {
    string initLine, word;
    vector<int> graphInfo;

    ifstream fileReader("graphs/" + fileName + ".txt");

    getline(fileReader, initLine);
    stringstream ss(initLine);

    while (ss >> word)
        graphInfo.push_back(stoi(word));

    vertexCount = graphInfo[0];
    edgeCount = graphInfo[1];

    vector<vector<int>> adj(vertexCount, vector<int>(vertexCount, INT_MAX));
    for (int i = 0; i < vertexCount; i++) adj[i][i] = 0;
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
