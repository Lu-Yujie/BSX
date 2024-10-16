#ifndef SUBGRAPHMATCHING_GENERATEQUERYPLAN_H
#define SUBGRAPHMATCHING_GENERATEQUERYPLAN_H

#include "graph/graph.h"
#include <vector>
#include <unordered_set>
class GenerateQueryPlan {
public:
    static void generateGQLQueryPlan(const Graph *data_graph, const Graph *query_graph, ui *candidates_count,
                                         ui *&order, ui *&pivot);

    static void checkQueryPlanCorrectness(const Graph* query_graph, ui* order, ui* pivot);

    static void printSimplifiedQueryPlan(const Graph* query_graph, ui* order);
private:
    static VertexID selectGQLStartVertex(const Graph *query_graph, ui *candidates_count);

    static void updateValidVertices(const Graph* query_graph, VertexID query_vertex, std::vector<bool>& visited,
                                    std::vector<bool>& adjacent);
};

#endif //SUBGRAPHMATCHING_GENERATEQUERYPLAN_H
