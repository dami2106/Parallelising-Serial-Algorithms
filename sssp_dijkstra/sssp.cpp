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

vector<vector<int> > makeGraph(int &vertexCount, int &edgeCount, const string &fileName);

void dijkstra(int vertexCount, vector<vector<int> > adj, vector<int> &dist, vector<int> &parents);

void checkInput(string graphName);

//The first argument argv[1] indicates which graph to run the algorithm on
int main(int argc, char *argv[]) {
    int vertexCount, edgeCount;    //Stores the variables needed for the algorithm
    double startTime, runTime = 0; //Needed to track the execution time

    checkInput(argv[1]); //Check the input graph is valid

    //Create the adjacency matrix and update the edge, vertex count variables for the given graph
    vector<vector<int> > adj = makeGraph(vertexCount, edgeCount, argv[1]);

    //2 vectors to store the distance array and parent array
    vector<int> dist(vertexCount, INT_MAX);
    vector<int> parents(vertexCount);

    //Run the SSSP algorithm and time it
    startTime = omp_get_wtime();
    dijkstra(vertexCount, adj, dist, parents);
    runTime = omp_get_wtime() - startTime;

    //Serial implementation has no speed up or validation so return the runtime
    cout << "Serial Time: " << runTime << endl;
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
void dijkstra(int vertexCount, vector<vector<int> > adj, vector<int> &l, vector<int> &parents) {
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
            if (adj[v][u] != INT_MAX)
                if (vT.find(v) == vT.end() && l[v] > l[u] + adj[v][u]) {
                    l[v] = l[u] + adj[v][u];
                    parents[v] = u;
                }
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

    for (int i = 0; i < vertexCount; i++)
        adj[i][i] = 0;

    //Close the file reader and return
    fileReader.close();
    return adj;
}
