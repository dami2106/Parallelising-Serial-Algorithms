#include <iostream>
#include <vector>
#include <sstream>
#include <array>
#include <unordered_set>
#include <omp.h>
#include <fstream>
#include <climits>

using namespace std;

vector<vector<int> > makeGraph(int &vertexCount, int &edgeCount, const string &fileName);
vector<int> dijkstra(int vertexCount, int startVertex, vector<vector<int> > adj);

void printPath(int vert, vector<int> parents);
void printSolution(int startVertex, vector<int> distances, vector<int> parents);

//arg[1] = file name, arg[2] = iterations
int main(int argc, char* argv[]) {
    int vertexCount, edgeCount, startVertex = 0;
    if(argv[1] == NULL) argv[1] = "graph_0"; //Fallback incase an argument is not given

    vector<vector<int> > adj = makeGraph(vertexCount, edgeCount, argv[1]);
    vector<int> distance = dijkstra(vertexCount, startVertex, adj);
    for (int c : distance)
        cout << c << " ";
}

vector<int> dijkstra(int vertexCount, int startVertex, vector<vector<int> > adj) {
    unordered_set<int> vT;
    vector<int> dist(vertexCount, INT_MAX);
    vector<int> parents(vertexCount);
    dist[startVertex] = 0;
    parents[startVertex] = -1;

    while (vT.size() != vertexCount) {
        int u, min = INT_MAX;

        for (int i = 0; i < vertexCount; i++)
            if (vT.find(i) == vT.end() && dist[i] < min) {
                min = dist[i];
                u = i;
            }

        vT.insert(u);

        for (int v = 0; v < vertexCount; v++) {
            if (adj[v][u] != -1)
                if (vT.find(v) == vT.end() && dist[v] > dist[u] + adj[v][u]) {
                    dist[v] = dist[u] + adj[v][u];
                    parents[v] = u;
                }
        }
    }
    return dist;
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
 * Creates and returns an adjacency matrix from std in
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

    vector<vector<int> > adj(vertexCount, vector<int>(vertexCount, -1));
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
