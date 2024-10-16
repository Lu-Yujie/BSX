#ifndef SUBGRAPHMATCHING_BUILDTABLE_H
#define SUBGRAPHMATCHING_BUILDTABLE_H

#include "graph/graph.h"
#include <vector>
class BuildEdgeIndex {
public:
    static void buildCansIdxIndex(const Graph* data_graph, const Graph* query_graph, ui** candidates, ui* candidates_count,
                            Edges*** edge_matrix);
    static void buildCansIndex(const Graph* data_graph, const Graph* query_graph, ui** candidates, ui* candidates_count,
                            Edges*** edge_matrix);
};

#endif //SUBGRAPHMATCHING_BUILDTABLE_H
