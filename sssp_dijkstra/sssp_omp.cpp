#include <iostream>
#include <vector>
#include <sstream>
#include <array>
#include <unordered_set>
#include <omp.h>
#include <fstream>
#include <climits>

using namespace std;

#define START 0
#define NUMTHREADS 8

vector<vector<int> > makeGraph(int &vertexCount, int &edgeCount, const string &fileName);

void serialDijkstra(int vertexCount, vector<vector<int> > adj, vector<int> &l, vector<int> &parents);

void parallelDijkstra(int vertexCount, vector<vector<int> > adj, vector<int> &l, vector<int> &parent);

void checkInput(string graphName);

//The first argument argv[1] indicates which graph to run the algorithm on
int main(int argc, char *argv[]) {
    int vertexCount, edgeCount; //Stores the variables needed for the algorithm (metadata about the current graph)
    double startTime, serRunTime = 0, parRunTime = 0; //Variables to keep track of the timing

    checkInput(argv[1]); //Check the input graph is valid

    //Create the adjacency matrix and update the edge, vertex count variables for the given graph
    vector<vector<int> > adj = makeGraph(vertexCount, edgeCount, argv[1]);

    //2 vectors to store the distance array and parent array of the serial code
    vector<int> serialDist(vertexCount, INT_MAX);
    vector<int> serialParents(vertexCount, -1);

    //2 vectors to store the distance array and parent array of the parallel code
    vector<int> parallelDist(vertexCount, INT_MAX);
    vector<int> parallelParents(vertexCount, -1);

    //Time and perform the serial SSSP dijkstra
    startTime = omp_get_wtime();
    serialDijkstra(vertexCount, adj, serialDist, serialParents);
    serRunTime += omp_get_wtime() - startTime;

    //Time and perform the parallel SSSP dijkstra
    startTime = omp_get_wtime();
    parallelDijkstra(vertexCount, adj, parallelDist, parallelParents);
    parRunTime += omp_get_wtime() - startTime;

    //Validate the parallel data against the serial distance array and serial parallel array
    if (serialDist != parallelDist) {
        cout << "(Validation Failed!)\n";
    } else {
        cout << "(Validation Passed!)\n";
        cout << "Serial Time : " << serRunTime << endl;
        cout << "Parallel Time : " << parRunTime << endl;
        cout << "Speed-Up : " << (serRunTime / parRunTime) << endl;
        cout << "Efficiency : " << (serRunTime / parRunTime) / NUMTHREADS << endl;
    }
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
    parents[START] = START; //Set the parent to the start vertex to the start vertex

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


/*
 * Implementation of OMP parallel SSSP Dijkstra. Performs the algorithm on the given graph while keeping track
 * of the min distances and path from START to all other nodes
 */
void parallelDijkstra(int vertexCount, vector<vector<int> > adj, vector<int> &l, vector<int> &parent) {
    //A vector representing the visited set to store the visited vertices
    vector<bool> vT(vertexCount, false);
    //Variables initialised for use in the parallel code.
    int threadID, threadCount, localMin = INT_MAX, localU = -1, currentVert, threadBoundLeft, threadBoundRight;
    int u = -1, min = INT_MAX;

    //Set the distance to the start vertex to 0
    l[START] = 0;
    parent[START] = START;

    //Create the parallel region while specifying the datascope of each of the above varaible
#pragma omp parallel num_threads(NUMTHREADS) firstprivate(localMin, localU) private(threadID, threadCount, currentVert, threadBoundLeft, threadBoundRight) shared(vT, min, u, adj, l, parent)
    {
        threadID = omp_get_thread_num(); //Stores the current thread number
        threadCount = omp_get_num_threads(); //Stores the number of threads

        //Sets the left & right bound of the current thread (used to ensure each thread only operates on its portion of the
        //distance array and adjacency list). Same effect as creating local arrays without the overhead (1D blocking)
        threadBoundLeft = threadID * (vertexCount / threadCount);
        threadBoundRight = ((threadID + 1) * (vertexCount / threadCount)) - 1;

        //Iterate for every vertex in the graph
        for (currentVert = 0; currentVert < vertexCount; currentVert++) {
//Allow a single thread to initialise and reset the global variables
#pragma omp single
            {
                u = -1;
                min = INT_MAX;
            }

            //Each thread resets its local closest vertex
            localU = -1;
            localMin = INT_MAX;

            //For each thread's subset of data, find the min/the closest vertex to the current vertex
            for (int i = threadBoundLeft; i <= threadBoundRight; i++) {
                if (!vT[i] && l[i] < localMin) localMin = l[i], localU = i;
            }

            //If the current threads min is closer than the global min, update the global min
            //1 thread at a time to avoid a race condition
#pragma omp critical
            {
                if (localMin < min) {
                    min = localMin;
                    u = localU;
                }
            }


#pragma omp barrier //Synchronise threads

#pragma omp single
            if (u != -1) vT[u] = true; //Set the current global closest vertex to visited by a single threads

#pragma omp barrier
            if (u != -1) {
                //For each vertex in the subset of vertices, update the current min distance and update the parent
                //array if the new distance is closer than the previous
                for (int i = threadBoundLeft; i <= threadBoundRight; i++) {
                    if (adj[i][u] != INT_MAX && l[u] != INT_MAX) {
                        if (!vT[i] && l[i] > l[u] + adj[i][u]) {
                            l[i] = l[u] + adj[i][u];
                            parent[i] = u;
                        }
                    }

                }
            }
#pragma omp barrier //Synchronise all threads again
        }
    }
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
