// compute min vertex cover
// not tested yet for survey graph format
#ifndef MIN_VERTEX_COVER_H
#define MIN_VERTEX_COVER_H

#include <vector>
#include <algorithm>
#include "graph/graph.h"
using namespace std;

class MinVCover {
public:
static ui minVertexCover(const Graph* graph, VertexID *order) {
    ui cover_num = 0;
    ui n = graph->getVerticesCount();
    auto offset = graph->getOffsets();
    auto edges = graph->getEdges();
    vector<int> degree(n, 0);
    vector<bool> visited(n, false);

    for (VertexID i = 0; i < n; i++) {
        if (graph->getVertexDegree(i) > 1)
            degree[i] = graph->getVertexDegree(i);
    }

    // select node with largest degree
    while (true) {
        int maxDegree = 0, maxNode = (VertexID)-1;
        for (VertexID i = 0; i < n; i++) {
            if (!visited[i] && degree[i] > maxDegree) {
                maxDegree = degree[i];
                maxNode = i;
            }
        }

        if (maxNode == (VertexID)-1) {
            break;  // all edges are covered
        }

        visited[maxNode] = true;
        order[cover_num++] = maxNode;

        // delete joint edges
        for (ui i = offset[maxNode]; i < offset[maxNode + 1]; i++) {
            VertexID neighbor = edges[i];
            degree[neighbor]--;
        }
    }

    return cover_num;
}
};

#endif