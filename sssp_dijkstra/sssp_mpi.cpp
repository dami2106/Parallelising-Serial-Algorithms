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

vector<int> parallelDijkstra(int vertexCount, int startVertex, vector<vector<int>> adj, vector<int> globalDist);

void printPath(int vert, vector<int> parents);

void printSolution(int startVertex, vector<int> distances, vector<int> parents);

int main(int argc, char *argv[]) {
    double startTime = 0, parRunTime = 0, serRunTime = 0;

    MPI_Init(NULL, NULL);

    int vertexCount, edgeCount, threadID, threadCount, flatSize, flatDist;
    vector<int> serDij;

    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCount);

//    vector<vector<int>> adj;
    vector<int> flatAdj;

    if (threadID == 0) {
        vector<vector<int>> adj = makeGraph(vertexCount, edgeCount, argv[1]);

        flatSize = vertexCount * vertexCount;
        flatDist = flatSize / threadCount;

        vector<int> temp(flatSize, INT_MAX);
        flatAdj = temp;

        for (int i = 0; i < vertexCount; i++) {
            for (int j = 0; j < vertexCount; j++) {
                flatAdj[i * vertexCount + j] = adj[i][j];
            }
        }

        startTime = MPI_Wtime();
        serDij = serialDijkstra(vertexCount, START, adj);
        serRunTime = MPI_Wtime() - startTime;
    }
    MPI_Bcast(&flatSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&flatDist, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&vertexCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<int> localAdj(flatDist, INT_MAX);
    MPI_Scatter(flatAdj.data(), flatDist, MPI_INT, localAdj.data(), flatDist, MPI_INT, 0, MPI_COMM_WORLD);

    vector<int> globalDist(vertexCount, 0);

    if (threadID == 0) startTime = MPI_Wtime();

    int localVals[2] = {INT_MAX, -1};
    int globalVals[2] = {INT_MAX, -1};

    int localCount = vertexCount / threadCount;
    int lowerBound = threadID * localCount;
    int upperBound = ((threadID + 1) * localCount) - 1;

    vector<bool> visited(localCount, false);
    vector<int> dist(localCount, INT_MAX);

    if (lowerBound <= START && START <= upperBound) {
        dist[START] = 0;
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
        //MPI_Barrier(MPI_COMM_WORLD);

        if (globalVals[1] == localVals[1]) {
            visited[(globalVals[1] - lowerBound)] = true;
        }


        for (int j = 0; j < localCount; j++) {
            int t = localAdj[j * vertexCount + globalVals[1]];
            if (t != INT_MAX) {
                if (!visited[j] && ((globalVals[0] + t) < dist[j])) {
                    dist[j] = globalVals[0] + t;
                }
            }

        }

    }

    if (threadID == 0)
        parRunTime = MPI_Wtime() - startTime;

    MPI_Gather(dist.data(), localCount, MPI_INT, globalDist.data(), localCount, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    if (threadID == 0) {
        if(globalDist != serDij)
            cout << "(Validation Failed!)";
        else
            cout << serRunTime/parRunTime << "  (Validation Passed!)";
    }
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
            if (adj[v][u] != INT_MAX)
                if (vT.find(v) == vT.end() && l[v] > l[u] + adj[v][u]) {
                    l[v] = l[u] + adj[v][u];
                    //parents[v] = u;
                }
        }
    }
    return l;
}

vector<int> parallelDijkstra(int vertexCount, int startVertex, vector<int> adj, vector<int> globalDist) {
//    int threadID, threadCount;
//    int localVals[2] = {-1, INT_MAX};
//    int globalVals[2] = {-1, INT_MAX};
//
//    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);
//    MPI_Comm_size(MPI_COMM_WORLD, &threadCount);
//
//    int localCount = vertexCount/threadCount;
//    int lowerBound = threadID * localCount;
//    int upperBound = ((threadID + 1) * localCount) - 1;

//    vector<bool> visited(localCount, false);
//    vector<int> dist(localCount, INT_MAX);
//
//    for(int i = 0 ; i < vertexCount ; i++)
//        dist[i] = adj[START * vertexCount + i];

//    if(lowerBound <= START && START <= upperBound)
//        visited[START - lowerBound] = true;
//
//    for(int i = 0 ; i < vertexCount ; i++) {
//        localVals[0] = -1;
//        localVals[1] = INT_MAX;
//
//        for(int j = 0 ; j < localCount ; j++){
//            if(visited[j] && dist[j] < localVals[1]){
//                localVals[0] = j;
//                localVals[1] = dist[j];
//            }
//        }
//
//        MPI_Allreduce(localVals, globalVals, 1, MPI_2INTEGER, MPI_MINLOC, MPI_COMM_WORLD);
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        if(globalVals[0] == localVals[0])
//            visited[localVals[0] - lowerBound] = true;
//
//        for(int j = 0 ; j < localCount ; j++){
//            if(visited[j] && (localVals[1] +  adj[j][localVals[0]]) < dist[j]){
//                dist[j] = (localVals[1] +  adj[j][localVals[0]]);
//            }
//        }
//        MPI_Gather(dist.data(), localCount, MPI_INT, globalDist.data(), localCount, MPI_INT, 0, MPI_COMM_WORLD);
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        if(threadID == 0) {
//            for (int c : globalDist) {
//                cout << c << " ";
//            }
//            cout << "\n";
//        }
    //return dist;
    //}
    vector<int> f;
    return f;
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
