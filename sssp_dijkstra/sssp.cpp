#include <iostream>
#include <vector>
#include <sstream>
#include <array>
#include <set>

using namespace std;

vector<vector<int> > makeGraph(int &vertexCount, int &edgeCount);

int main() {
    int vertexCount, edgeCount;
    vector<vector<int> > adj = makeGraph(vertexCount, edgeCount);;
    set<int> visited;
    for(int i = 1 ; i < vertexCount ; i++) {
        
    }


}

/*
 * Creates and returns an adjacency matrix from std in
 */
vector<vector<int> > makeGraph(int &vertexCount, int &edgeCount) {
    string initLine, word;
    vector<int> graphInfo;

    getline(cin, initLine);
    stringstream ss(initLine);

    while(ss >> word)
        graphInfo.push_back(stoi(word));

    vertexCount = graphInfo[0];
    edgeCount = graphInfo[1];

    vector< vector<int> > adj(vertexCount, vector<int>(vertexCount, -1));
    for (int i = 0; i < edgeCount; ++i) {
        string currNode, nodeItem;
        array<int, 3> nodeInfo;
        int k = 0;

        getline(cin, currNode);
        stringstream sx(currNode);

        while (sx >> nodeItem){
            nodeInfo[k] = stoi(nodeItem);
            ++k;
        }
        adj[nodeInfo[0]][nodeInfo[1]] = nodeInfo[2];
    }

    return adj;
}
