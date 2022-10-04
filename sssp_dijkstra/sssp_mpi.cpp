#include <iostream>
#include <vector>
#include <sstream>
#include <array>
#include <unordered_set>
#include <fstream>
#include <climits>
#include <mpi.h>

using namespace std;

#define START 0

vector<vector<int>> makeGraph(int &vertexCount, int &edgeCount, const string &fileName);

vector<int> serialDijkstra(int vertexCount, int startVertex, vector<vector<int>> adj);

vector<int> parallelDijkstra(int vertexCount, int startVertex, vector<vector<int>> adj);

void printPath(int vert, vector<int> parents);

void printSolution(int startVertex, vector<int> distances, vector<int> parents);

int main(int argc, char *argv[]) {
    vector<vector<int>> adj;
    double startTime, serRunTime = 0, parRunTime = 0;
    int vertexCount, edgeCount, startVertex = 0, id;
    int iterations = atoi(argv[2]);

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if(id == 0)
        adj = makeGraph(vertexCount, edgeCount, argv[1]);

    MPI_Bcast(adj.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    vector<int> dist = parallelDijkstra(vertexCount, startVertex, adj);

    cout << "Serial Time : " << serRunTime / iterations << endl;
    cout << "Parallel Time : " << parRunTime / iterations << endl;
    cout << "Speed-Up : " << (serRunTime / parRunTime) << endl;
}

vector<int> serialDijkstra(int vertexCount, int startVertex, vector<vector<int>> adj) {
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
            if (adj[v][u] != -1)
                if (vT.find(v) == vT.end() && l[v] > l[u] + adj[v][u]) {
                    l[v] = l[u] + adj[v][u];
                    //parents[v] = u;
                }
        }
    }
    return l;
}

vector<int> parallelDijkstra(int vertexCount, int startVertex, vector<vector<int>> adj) {
    int threadID, threadCount;
    int localVals[2] = {-1, INT_MAX};
    int globalVals[2] = {-1, INT_MAX};

    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCount);

    int localCount = vertexCount/threadCount;
    int lowerBound = threadID * localCount;
    int upperBound = ((threadID + 1) * localCount) - 1;

    vector<bool> visited(localCount, false);
    vector<int> dist(localCount, INT_MAX);

    if(threadID == 0) visited[0] = true;

    for(int i = 0 ; i < vertexCount ; i++) {
        localVals[0] = -1;
        localVals[1] = INT_MAX;

        for(int j = 0 ; j < localCount ; j++){
            if(visited[j] && dist[j] < )
        }

    }
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

    vector<vector<int>> adj(vertexCount, vector<int>(vertexCount, -1));
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
