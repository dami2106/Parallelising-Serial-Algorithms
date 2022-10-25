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

void serialDijkstra(int vertexCount, vector<vector<int> > adj, vector<int> &l, vector<int> &parents);

void parallelDijkstra(int vertexCount, vector<int> localAdj, vector<int> &globalDist, vector<int> &globalParent);

void checkInput(string graphName);


//The first argument argv[1] indicates which graph to run the algorithm on
int main(int argc, char *argv[]) {
    //Variables to keep track of the timing
    double startTime = 0, parRunTime = 0, serRunTime = 0;
    //Stores the variables needed for the algorithm (metadata about the current graph)
    //Also stores varaibles needed for MPI to prepare the data
    int vertexCount, edgeCount, threadID, threadCount, flatSize, flatDist;

    //Initialise MPI
    MPI_Init(NULL, NULL);
    //Get the current threads rank and the number of threads
    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCount);

    vector<int> flatAdj; //A vector that stores a flattened adjacency matrix, so it can be scattered
    vector<vector<int>> adj; //Stores the adjacency matrix

    //Use a single thread to do the initialisation
    if (threadID == 0) {

        checkInput(argv[1]); //Check the input graph is valid

        //Setup and create the adjacency matrix
        adj = makeGraph(vertexCount, edgeCount, argv[1]);

        //Determine the size of the flatened adjacency matrix
        flatSize = vertexCount * vertexCount;
        //Determine the size of the distributed (scattered) adjacency matrix
        flatDist = flatSize / threadCount;

        //A temporary vertex used to initialise flatAdj
        vector<int> temp(flatSize, INT_MAX);
        flatAdj = temp;

        //For each vertex in the adjacency list, save its value to the flattened list
        for (int i = 0; i < vertexCount; i++) {
            for (int j = 0; j < vertexCount; j++) {
                flatAdj[i * vertexCount + j] = adj[i][j];
            }
        }
    }

    //Broadcast the size of the flat adjacency list, the size of the distributed list and the number of vertices
    //From thread 0 to all the other threads
    MPI_Bcast(&flatSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&flatDist, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&vertexCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //Arrays needed to store the local adjacency list for each thread
    //as well as the output arrays of the parallel implementation - distance and parents
    vector<int> localAdj(flatDist, INT_MAX);
    vector<int> parallelDist(vertexCount, INT_MAX);
    vector<int> parallelParent(vertexCount, 0);

    //Arrays needed to store the output of the serial implementation
    vector<int> serialDist(vertexCount, INT_MAX);
    vector<int> serialParent(vertexCount, 0);

    //Scatter the flattened adjacency matrix to each thread before running the parallel code
    //Scatters before the code is called to ensure fair computation
    MPI_Scatter(flatAdj.data(), flatDist, MPI_INT, localAdj.data(), flatDist, MPI_INT, 0, MPI_COMM_WORLD);

    //Synchronise all threads
    MPI_Barrier(MPI_COMM_WORLD);

    //Use thread 0 time get the time and validation data from the serial implementation
    if (threadID == 0) {
        startTime = MPI_Wtime();
        serialDijkstra(vertexCount, adj, serialDist, serialParent);
        serRunTime = MPI_Wtime() - startTime;
    }

    //Synchronise all threads
    MPI_Barrier(MPI_COMM_WORLD);

    //Use thread 0 to time the parallel implementation as well as running the parallel implementation
    if (threadID == 0) startTime = MPI_Wtime();
    parallelDijkstra(vertexCount, localAdj, parallelDist, parallelParent);
    if (threadID == 0) parRunTime = MPI_Wtime() - startTime;

    //Synchronise all threads
    MPI_Barrier(MPI_COMM_WORLD);


    if (threadID == 0) {
        if ((serialDist != parallelDist) || (serialParent != parallelParent))
            cout << "(Validation Failed!)";
        else {
            cout << "(Validation Passed!)\n";
            cout << "Serial Time : " << serRunTime << endl;
            cout << "Parallel Time : " << parRunTime << endl;
            cout << "Speed-Up : " << (serRunTime / parRunTime) << endl;
            cout << "Efficiency : " << (serRunTime / parRunTime) / threadCount << endl;
        }
    }
    MPI_Finalize();
}

/*
 * Checks if the provided argument is correct, exits if it isn't
 */
void checkInput(string graphName) {
    if (graphName.length() != 7) cout << "INCORRECT GRAPH NAME\n", exit(1);
    if (graphName[5] != '_') cout << "INCORRECT GRAPH NAME\n", exit(1);
    string graphNum(1, graphName[6]);
    if (stoi(graphNum) < 0 || stoi(graphNum) > 7) cout << "INCORRECT GRAPH NAME\n", exit(1);
}

/*
 * Implementation of serial SSSP Dijkstra. Performs the algorithm on the given graph while keeping track
 * of the min distances and path from START to all other nodes
 */
void serialDijkstra(int vertexCount, vector<vector<int> > adj, vector<int> &l, vector<int> &parents) {
    unordered_set<int> vT; //Set to store visited vertices
    l[START] = 0; //Set the distance to the start vertex to 0

    while ((int) vT.size() != vertexCount) {
        int u = -1, min = INT_MAX;

        //Iterate for every vertex
        for (int i = 0; i < vertexCount; i++)
            //Find the closest vertex with the smallest distance that hasn't been visited
            if (vT.find(i) == vT.end() && l[i] < min) {
                min = l[i];
                u = i;
            }

        //Add the current closest vertex to the visited list
        vT.insert(u);

        //For every vertex again
        for (int v = 0; v < vertexCount; v++) {
            //Check if the new path to this vertex is shorter than the previous one stored
            //If it is, update the distance and set the new parent of this vertex
            if (adj[v][u] != INT_MAX) {
                if (vT.find(v) == vT.end() && l[v] > l[u] + adj[v][u]) {
                    l[v] = l[u] + adj[v][u];
                    parents[v] = u;
                }
            }
        }
    }
}

void parallelDijkstra(int vertexCount, vector<int> localAdj, vector<int> &globalDist, vector<int> &globalParent) {
    int threadCount, threadID; //Store the current thread rank and thread ID
    //Use an array to store the local min vertex and min distance (so reduction can be done in 1 step)
    int localVals[2] = {INT_MAX, -1};
    //Use an array to store the global min vertex and min distance (so reduction can be done in 1 step)
    int globalVals[2] = {INT_MAX, -1};

    //Get the thread rank and ID
    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCount);

    //Get the number of vertices given to the current thread as well as the lower bound and upper bound
    //which will be used to index the flattened adjacency matrix
    int localCount = vertexCount / threadCount;
    int lowerBound = threadID * localCount;
    int upperBound = ((threadID + 1) * localCount) - 1;


    vector<bool> visited(localCount, false); //Have a vector that stores which local vertices have been visited
    vector<int> dist(localCount, INT_MAX); //Vector that stores the distances for the local vertices
    vector<int> parent(localCount, -1);//Vector that stores the parents for the local vertices

    //Set the distance to the start vertex equal to 0
    if (lowerBound <= START && START <= upperBound) {
        dist[START] = 0;
        parent[START] = 0;
    }

    //For each vertex in the graph
    for (int i = 0; i < vertexCount; i++) {
        localVals[0] = INT_MAX; //Reset the local min vertex and min distance
        localVals[1] = -1;

        //For each local vertex in the thread
        for (int j = 0; j < localCount; j++) {
            //Find the min distance and closest vertex to the current vertex
            if (!visited[j] && (dist[j] < localVals[0])) {
                localVals[1] = j + lowerBound;
                localVals[0] = dist[j];
            }
        }
        //Reduce the closest vertex to all threads after finding the min
        MPI_Allreduce(localVals, globalVals, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);

        //If the current thread holds the closest vertex
        if (globalVals[1] == localVals[1]) {
            visited[(globalVals[1] - lowerBound)] = true; //Set the thread to be visited within this threaD
        }

        //For each local vertex in the current thread
        for (int j = 0; j < localCount; j++) {
            //Compute the new possible distance
            int t = localAdj[j * vertexCount + globalVals[1]];
            if (t != INT_MAX) {
                //Check if the new distance is closer and if it is update the
                // //current distance to this vertex and the parent
                if (!visited[j] && ((globalVals[0] + t) < dist[j])) {
                    dist[j] = globalVals[0] + t;
                    parent[j] = globalVals[1];
                }
            }

        }

    }

    //Gather the local distance and parent arrays to the global arrays into thread 0
    MPI_Gather(dist.data(), localCount, MPI_INT, globalDist.data(), localCount, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(parent.data(), localCount, MPI_INT, globalParent.data(), localCount, MPI_INT, 0, MPI_COMM_WORLD);
}


/*
 * Creates and returns an adjacency matrix from std in
 */
vector<vector<int> > makeGraph(int &vertexCount, int &edgeCount, const string &fileName) {
    string initLine, word;
    vector<int> graphInfo;

    //Read in the graph from a text file with ifstream
    ifstream fileReader("graphs/" + fileName + ".txt");

    //Read the file into a string and split it up using string stream
    getline(fileReader, initLine);
    stringstream ss(initLine);

    while (ss >> word)
        graphInfo.push_back(stoi(word));

    //Set up the vertex and edge count
    vertexCount = graphInfo[0];
    edgeCount = graphInfo[1];

    //Initialise a 2D adjacency matrix with INT_MAX everywhere
    vector<vector<int> > adj(vertexCount, vector<int>(vertexCount, INT_MAX));
    //Iterate for each edge
    for (int i = 0; i < edgeCount; ++i) {
        string currNode, nodeItem;
        array<int, 3> nodeInfo;
        int k = 0;

        //Read in the edge vertices and the weight
        getline(fileReader, currNode);
        stringstream sx(currNode);

        while (sx >> nodeItem) {
            nodeInfo[k] = stoi(nodeItem);
            ++k;
        }

        //Set the weight between the 2 vertices in the adjacency matrix & also make sure it stays symmetric
        adj[nodeInfo[0]][nodeInfo[1]] = nodeInfo[2];
        adj[nodeInfo[1]][nodeInfo[0]] = nodeInfo[2];
    }

    //Close the file reader and return
    fileReader.close();
    return adj;
}