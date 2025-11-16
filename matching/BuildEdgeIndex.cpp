#include "BuildEdgeIndex.h"
#include <vector>
#include <algorithm>
using namespace std;

void BuildEdgeIndex::buildCansIdxIndex(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count,
                             Edges ***edge_matrix) {
    auto q_num = query_graph->getVerticesCount();
    auto d_num = data_graph->getVerticesCount();
    ui* flag = new ui[d_num];
    ui* updated_flag = new ui[d_num];
    fill(flag, flag + d_num, 0);

    for (ui i = 0; i < q_num; ++i) {
        for (ui j = 0; j < q_num; ++j) {
            edge_matrix[i][j] = nullptr;
        }
    }

    // generate table order based on node degree
    vector<VertexID> build_table_order(q_num);
    for (ui i = 0; i < q_num; ++i) {
        build_table_order[i] = i;
    }

    sort(build_table_order.begin(), build_table_order.end(), [query_graph](VertexID l, VertexID r) {
        if (query_graph->getVertexDegree(l) == query_graph->getVertexDegree(r)) {
            return l < r;
        }
        return query_graph->getVertexDegree(l) > query_graph->getVertexDegree(r);
    });

    vector<ui> temp_edges(data_graph->getEdgesCount() * 2);

    for (auto u : build_table_order) {
        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
#ifdef ELABELED_GRAPH
        auto u_elabels = query_graph->getVertexEdgeLabels(u, u_nbrs_count);
#endif
        ui updated_flag_count = 0;

        for (ui i = 0; i < u_nbrs_count; ++i) {
            VertexID u_nbr = u_nbrs[i];
#ifdef ELABELED_GRAPH
            LabelID u_elabel = u_elabels[i];
#endif
            if (edge_matrix[u][u_nbr] != nullptr)
                continue;

            // u--node on q_graph  v--candidates
            // updated_flag_count: #candidates or #updated_flag has been set
            // flag: idx of candidate, flag[v]!=0 means v is candidate of u
            // updated_flag:  recover flag, because reconstruct one is too large
            if (updated_flag_count == 0) {
                for (ui j = 0; j < candidates_count[u]; ++j) {
                    VertexID v = candidates[u][j];
                    flag[v] = j + 1;
                    updated_flag[updated_flag_count++] = v;
                }
            }

            edge_matrix[u_nbr][u] = new Edges;
            edge_matrix[u_nbr][u]->vertex_count_ = candidates_count[u_nbr];
            edge_matrix[u_nbr][u]->offset_ = new ui[candidates_count[u_nbr] + 1];

            edge_matrix[u][u_nbr] = new Edges;
            edge_matrix[u][u_nbr]->vertex_count_ = candidates_count[u];
            edge_matrix[u][u_nbr]->offset_ = new ui[candidates_count[u] + 1];
            fill(edge_matrix[u][u_nbr]->offset_, edge_matrix[u][u_nbr]->offset_ + candidates_count[u] + 1, 0);

            ui local_edge_count = 0;
            ui local_max_degree = 0;

            for (ui j = 0; j < candidates_count[u_nbr]; ++j) {
                VertexID v = candidates[u_nbr][j];
                edge_matrix[u_nbr][u]->offset_[j] = local_edge_count;

                ui v_nbrs_count;
                const VertexID* v_nbrs = data_graph->getVertexNeighbors(v, v_nbrs_count);
#ifdef ELABELED_GRAPH
                auto v_elabels = data_graph->getVertexEdgeLabels(v, v_nbrs_count);
#endif
                ui local_degree = 0;

                for (ui k = 0; k < v_nbrs_count; ++k) {
                    VertexID v_nbr = v_nbrs[k];
#ifdef ELABELED_GRAPH
                    LabelID v_elabel = v_elabels[k];
                    if (v_elabel != u_elabel) {
                        continue;
                    }
#endif
                    if (flag[v_nbr] != 0) {
                        ui position = flag[v_nbr] - 1;
                        temp_edges[local_edge_count++] = position;
                        edge_matrix[u][u_nbr]->offset_[position + 1] += 1;
                        local_degree += 1;
                    }
                }

                if (local_degree > local_max_degree) {
                    local_max_degree = local_degree;
                }
            }

            edge_matrix[u_nbr][u]->offset_[candidates_count[u_nbr]] = local_edge_count;
            edge_matrix[u_nbr][u]->max_degree_ = local_max_degree;
            edge_matrix[u_nbr][u]->edge_count_ = local_edge_count;
            edge_matrix[u_nbr][u]->edge_ = new ui[local_edge_count];
            copy(temp_edges.begin(), temp_edges.begin() + local_edge_count, edge_matrix[u_nbr][u]->edge_);

            edge_matrix[u][u_nbr]->edge_count_ = local_edge_count;
            edge_matrix[u][u_nbr]->edge_ = new ui[local_edge_count];

            local_max_degree = 0;
            for (ui j = 1; j <= candidates_count[u]; ++j) {
                if (edge_matrix[u][u_nbr]->offset_[j] > local_max_degree) {
                    local_max_degree = edge_matrix[u][u_nbr]->offset_[j];
                }
                edge_matrix[u][u_nbr]->offset_[j] += edge_matrix[u][u_nbr]->offset_[j - 1];
            }

            edge_matrix[u][u_nbr]->max_degree_ = local_max_degree;

            for (ui j = 0; j < candidates_count[u_nbr]; ++j) {
                ui begin = j;
                for (ui k = edge_matrix[u_nbr][u]->offset_[begin]; k < edge_matrix[u_nbr][u]->offset_[begin + 1]; ++k) {
                    ui end = edge_matrix[u_nbr][u]->edge_[k];

                    edge_matrix[u][u_nbr]->edge_[edge_matrix[u][u_nbr]->offset_[end]++] = begin;
                }
            }

            for (ui j = candidates_count[u]; j >= 1; --j) {
                edge_matrix[u][u_nbr]->offset_[j] = edge_matrix[u][u_nbr]->offset_[j - 1];
            }
            edge_matrix[u][u_nbr]->offset_[0] = 0;
        }

        for (ui i = 0; i < updated_flag_count; ++i) {
            VertexID v = updated_flag[i];
            flag[v] = 0;
        }
    }
}

void
BuildEdgeIndex::buildCansIndex(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count,
                               Edges ***edge_matrix) {
    auto q_num = query_graph->getVerticesCount();
    auto d_num = data_graph->getVerticesCount();
    ui* flag = new ui[d_num];
    ui* updated_flag = new ui[d_num];
    fill(flag, flag + d_num, 0);

    for (ui i = 0; i < q_num; ++i) {
        for (ui j = 0; j < q_num; ++j) {
            edge_matrix[i][j] = nullptr;
        }
    }

    vector<VertexID> build_table_order(q_num);
    for (ui i = 0; i < q_num; ++i) {
        build_table_order[i] = i;
    }

    sort(build_table_order.begin(), build_table_order.end(), [query_graph](VertexID l, VertexID r) {
        if (query_graph->getVertexDegree(l) == query_graph->getVertexDegree(r)) {
            return l < r;
        }
        return query_graph->getVertexDegree(l) > query_graph->getVertexDegree(r);
    });

    vector<ui> temp_edges(data_graph->getEdgesCount() * 2);

    for (auto u : build_table_order) {
        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
#ifdef ELABELED_GRAPH
        auto u_elabels = query_graph->getVertexEdgeLabels(u, u_nbrs_count);
#endif
        ui updated_flag_count = 0;

        for (ui i = 0; i < u_nbrs_count; ++i) {
            VertexID u_nbr = u_nbrs[i];
#ifdef ELABELED_GRAPH
            LabelID u_elabel = u_elabels[i];
#endif
            if (edge_matrix[u][u_nbr] != nullptr)
                continue;

            if (updated_flag_count == 0) {
                for (ui j = 0; j < candidates_count[u]; ++j) {
                    VertexID v = candidates[u][j];
                    flag[v] = j + 1;
                    updated_flag[updated_flag_count++] = v;
                }
            }

            edge_matrix[u_nbr][u] = new Edges;
            edge_matrix[u_nbr][u]->vertex_count_ = candidates_count[u_nbr];
            edge_matrix[u_nbr][u]->offset_ = new ui[candidates_count[u_nbr] + 1];

            edge_matrix[u][u_nbr] = new Edges;
            edge_matrix[u][u_nbr]->vertex_count_ = candidates_count[u];
            edge_matrix[u][u_nbr]->offset_ = new ui[candidates_count[u] + 1];
            fill(edge_matrix[u][u_nbr]->offset_, edge_matrix[u][u_nbr]->offset_ + candidates_count[u] + 1, 0);

            ui local_edge_count = 0;
            ui local_max_degree = 0;

            for (ui j = 0; j < candidates_count[u_nbr]; ++j) {
                VertexID v = candidates[u_nbr][j];
                edge_matrix[u_nbr][u]->offset_[j] = local_edge_count;

                ui v_nbrs_count;
                const VertexID* v_nbrs = data_graph->getVertexNeighbors(v, v_nbrs_count);
#ifdef ELABELED_GRAPH
                auto v_elabels = data_graph->getVertexEdgeLabels(v, v_nbrs_count);
#endif
                ui local_degree = 0;

                for (ui k = 0; k < v_nbrs_count; ++k) {
                    VertexID v_nbr = v_nbrs[k];
#ifdef ELABELED_GRAPH
                    LabelID v_elabel = v_elabels[k];
                    if (v_elabel != u_elabel) {
                        continue;
                    }
#endif
                    if (flag[v_nbr] != 0) {
                        ui position = flag[v_nbr] - 1;
                        temp_edges[local_edge_count++] = v_nbr;  // record id. instead of idx
                        edge_matrix[u][u_nbr]->offset_[position + 1] += 1;
                        local_degree += 1;
                    }
                }

                if (local_degree > local_max_degree) {
                    local_max_degree = local_degree;
                }
            }

            edge_matrix[u_nbr][u]->offset_[candidates_count[u_nbr]] = local_edge_count;
            edge_matrix[u_nbr][u]->max_degree_ = local_max_degree;
            edge_matrix[u_nbr][u]->edge_count_ = local_edge_count;
            edge_matrix[u_nbr][u]->edge_ = new ui[local_edge_count];
            copy(temp_edges.begin(), temp_edges.begin() + local_edge_count, edge_matrix[u_nbr][u]->edge_);

            edge_matrix[u][u_nbr]->edge_count_ = local_edge_count;
            edge_matrix[u][u_nbr]->edge_ = new ui[local_edge_count];

            local_max_degree = 0;
            for (ui j = 1; j <= candidates_count[u]; ++j) {
                if (edge_matrix[u][u_nbr]->offset_[j] > local_max_degree) {
                    local_max_degree = edge_matrix[u][u_nbr]->offset_[j];
                }
                edge_matrix[u][u_nbr]->offset_[j] += edge_matrix[u][u_nbr]->offset_[j - 1];
            }

            edge_matrix[u][u_nbr]->max_degree_ = local_max_degree;

            for (ui j = 0; j < candidates_count[u_nbr]; ++j) {
                VertexID begin = candidates[u_nbr][j];
                for (ui k = edge_matrix[u_nbr][u]->offset_[j]; k < edge_matrix[u_nbr][u]->offset_[j + 1]; ++k) {
                    ui end_idx = flag[edge_matrix[u_nbr][u]->edge_[k]] - 1;
                    edge_matrix[u][u_nbr]->edge_[edge_matrix[u][u_nbr]->offset_[end_idx]++] = begin;
                }
            }

            for (ui j = candidates_count[u]; j >= 1; --j) {
                edge_matrix[u][u_nbr]->offset_[j] = edge_matrix[u][u_nbr]->offset_[j - 1];
            }
            edge_matrix[u][u_nbr]->offset_[0] = 0;
        }

        for (ui i = 0; i < updated_flag_count; ++i) {
            VertexID v = updated_flag[i];
            flag[v] = 0;
        }
    }
}
