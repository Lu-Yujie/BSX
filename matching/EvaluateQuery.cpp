#include "EvaluateQuery.h"
#include "computesetintersection.h"
#include "rapidMatch/execution_tree/execution_tree_generator.h"
#include "bsx/IndepSet.h"
#include "bsx/nodeSim.h"
#include "bsx/SetOp.h"
#include "bsx/vcover.h"
#include "timeOp.h"
#include <vector>
#include <cstring>

#include "pretty_print.h"

std::function<bool(std::pair<std::pair<VertexID, ui>, ui>, std::pair<std::pair<VertexID, ui>, ui>)>
    EvaluateQuery::extendable_vertex_compare = [](std::pair<std::pair<VertexID, ui>, ui> l, std::pair<std::pair<VertexID, ui>, ui> r) {
    if (l.first.second == 1 && r.first.second != 1) {
        return true;
    }
    else if (l.first.second != 1 && r.first.second == 1) {
        return false;
    }
    else
    {
        return l.second > r.second;
    }
};

// 按照给定的order计算每个点邻居中匹配顺序在其之前的点的额数量
// 其中需要把pivot点去掉，（不懂）
void EvaluateQuery::generateBN(const Graph *query_graph, ui *order, ui *pivot, ui **&bn, ui *&bn_count) {
    ui q_num = query_graph->getVerticesCount();
    bn_count = new ui[q_num];
    std::fill(bn_count, bn_count + q_num, 0);
    bn = new ui *[q_num];
    for (ui i = 0; i < q_num; ++i) {
        bn[i] = new ui[q_num];
    }

    std::vector<bool> visited_vertices(q_num, false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < q_num; ++i) {
        VertexID vertex = order[i];

        ui nbrs_cnt;
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (ui j = 0; j < nbrs_cnt; ++j) {
            VertexID nbr = nbrs[j];

            if (visited_vertices[nbr] && nbr != pivot[i]) {
                bn[i][bn_count[i]++] = nbr;
            }
        }

        visited_vertices[vertex] = true;
    }
}

// 这个函数不需要去掉pivot（不懂）
void EvaluateQuery::generateBN(const Graph *query_graph, ui *order, ui **&bn, ui *&bn_count) {
    ui q_num = query_graph->getVerticesCount();
    bn_count = new ui[q_num];
    std::fill(bn_count, bn_count + q_num, 0);
    bn = new ui *[q_num];
    for (ui i = 0; i < q_num; ++i) {
        bn[i] = new ui[q_num];
    }

    std::vector<bool> visited_vertices(q_num, false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < q_num; ++i) {
        VertexID vertex = order[i];

        ui nbrs_cnt;
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (ui j = 0; j < nbrs_cnt; ++j) {
            VertexID nbr = nbrs[j];

            if (visited_vertices[nbr]) {
                bn[i][bn_count[i]++] = nbr;
            }
        }

        visited_vertices[vertex] = true;
    }
}

bool
EvaluateQuery::ExploreEngine(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                            ui *candidates_count, ui *order, ui *pivot, uint64_t output_limit_num, uint64_t &call_cnt,
                            mpz_t embedding_cnt, int64_t& time_limit) {
    // Generate the bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, pivot, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    // Evaluate the query.
    bool overtime = false;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidate_idx[cur_depth][i] = i;
    }

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            if (TimeOp::getClockNan() >= time_limit) {
                overtime = true;
                goto EXIT;
            }
            if ((ui)cur_depth == max_depth - 1) {
                mpz_add_ui(embedding_cnt, embedding_cnt, 1);
                visited_vertices[v] = false;   
                if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
                    goto EXIT;
                }
            } else {
                call_cnt += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidateIndex(data_graph, cur_depth, embedding, idx_embedding, idx_count,
                                            valid_candidate_idx, edge_matrix, visited_vertices, bn,
                                            bn_count, order, pivot, candidates, query_graph);
            }
        }

        // 回溯部分
        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
    }


    // Release the buffer.
    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    return overtime;
}

void
EvaluateQuery::allocateBuffer(const Graph *data_graph, const Graph *query_graph, ui *candidates_count, ui *&idx,
                              ui *&idx_count, ui *&embedding, ui *&idx_embedding, ui *&temp_buffer,
                              ui **&valid_candidate_idx, bool *&visited_vertices) {
    ui q_num = query_graph->getVerticesCount();
    ui d_num = data_graph->getVerticesCount();
    ui max_candidates_num = candidates_count[0];

    for (ui i = 1; i < q_num; ++i) {
        VertexID cur_vertex = i;
        ui cur_candidate_num = candidates_count[cur_vertex];

        if (cur_candidate_num > max_candidates_num) {
            max_candidates_num = cur_candidate_num;
        }
    }

    idx = new ui[q_num];
    idx_count = new ui[q_num];
    embedding = new ui[q_num];
    idx_embedding = new ui[q_num];
    visited_vertices = new bool[d_num];
    temp_buffer = new ui[max_candidates_num];
    valid_candidate_idx = new ui *[q_num];
    for (ui i = 0; i < q_num; ++i) {
        valid_candidate_idx[i] = new ui[max_candidates_num];
    }

    std::fill(visited_vertices, visited_vertices + d_num, false);
}

// 检查indices结构上相邻candidate的每一条边，返回所有边都能匹配上的结果，也就是下一个匹配点
void EvaluateQuery::generateValidCandidateIndex(const Graph *data_graph, ui depth, ui *embedding, ui *idx_embedding,
                                                ui *idx_count, ui **valid_candidate_index, Edges ***edge_matrix,
                                                bool *visited_vertices, ui **bn, ui *bn_cnt, ui *order, ui *pivot,
                                                ui **candidates, const Graph *query_graph) {
    VertexID u = order[depth];
    VertexID pivot_vertex = pivot[depth];
    ui idx_id = idx_embedding[pivot_vertex];
    Edges &edge = *edge_matrix[pivot_vertex][u];
    ui count = edge.offset_[idx_id + 1] - edge.offset_[idx_id];
    ui *candidate_idx = edge.edge_ + edge.offset_[idx_id];

    ui valid_candidate_index_count = 0;

    if (bn_cnt[depth] == 0) {
        for (ui i = 0; i < count; ++i) {
            ui temp_idx = candidate_idx[i];
            VertexID temp_v = candidates[u][temp_idx];

            if (!visited_vertices[temp_v])
                valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
        }
    } else {
        for (ui i = 0; i < count; ++i) {
            ui temp_idx = candidate_idx[i];
            VertexID temp_v = candidates[u][temp_idx];

            if (!visited_vertices[temp_v]) {
                bool valid = true;

                for (ui j = 0; j < bn_cnt[depth]; ++j) {
                    VertexID u_bn = bn[depth][j];
                    VertexID u_bn_v = embedding[u_bn];
#ifdef ELABELED_GRAPH
                    LabelID elabel = query_graph->getEdgeLabel(u, u_bn, true);
                    if (!data_graph->checkEdgeExistence(temp_v, u_bn_v, elabel)) {
#else
                    if (!data_graph->checkEdgeExistence(temp_v, u_bn_v)) {
#endif
                        valid = false;
                        break;
                    }
                }

                if (valid)
                    valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
            }
        }
    }

    idx_count[depth] = valid_candidate_index_count;
}

void EvaluateQuery::releaseBuffer(ui q_num, ui *idx, ui *idx_count, ui *embedding, ui *idx_embedding,
                                  ui *temp_buffer, ui **valid_candidate_idx, bool *visited_vertices, ui **bn,
                                  ui *bn_count) {
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] idx_embedding;
    delete[] visited_vertices;
    delete[] bn_count;
    delete[] temp_buffer;
    for (ui i = 0; i < q_num; ++i) {
        delete[] valid_candidate_idx[i];
        delete[] bn[i];
    }

    delete[] valid_candidate_idx;
    delete[] bn;
}

bool
EvaluateQuery::LFTJ(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                    ui *candidates_count, ui *order, uint64_t output_limit_num, uint64_t &call_cnt, uint64_t &valid_vtx_cnt,
                    mpz_t embedding_cnt, int64_t& time_limit) {
    // Generate bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    bool overtime = false;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    // // Valid vertex count statistics
    // bool ** valid_vertices = new bool* [max_depth];
    // for (int i = 0; i < max_depth; ++i) {
    //     valid_vertices[i] = new bool[data_vertices_count];
    //     std::fill(valid_vertices[i], valid_vertices[i] + data_vertices_count, false);
    // }

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidate_idx[cur_depth][i] = i;
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif
    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];

            if (visited_vertices[v]) {
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif
            if (TimeOp::getClockNan() >= time_limit) {
                overtime = true;
                goto EXIT;
            }
            if ((ui)cur_depth == max_depth - 1) {
                // for (ui d = 0; d < max_depth; ++d){
                //     valid_vertices[order[d]][embedding[d]] = true;
                // }
                mpz_add_ui(embedding_cnt, embedding_cnt, 1);
                visited_vertices[v] = false;

#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
                    goto EXIT;
                }
            } else {
                call_cnt += 1;
                cur_depth += 1;

                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);

#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0) {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                } else {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }
        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0) {
                if (!vec_failing_set[cur_depth].test(u)) {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                } else {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;
        }
    }

    // // Counting the valid vertices
    // for (ui i = 0; i < data_vertices_count; ++i){
    //     for (ui j = 0; j < max_depth; ++j){
    //         if (valid_vertices[j][i]){
    //             ++temp_valid_vtx_cnt;
    //         }
    //     }
    // }
    // valid_vtx_cnt = temp_valid_vtx_cnt;

    // Release the buffer.

    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    // // release the counting buffer
    // for (ui i = 0; i < max_depth; ++i) {;
    //     delete[] valid_vertices[i];
    // }
    // delete[] valid_vertices;

    return overtime;
}

void EvaluateQuery::generateValidCandidateIndex(ui depth, ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                                Edges ***edge_matrix, ui **bn, ui *bn_cnt, ui *order,
                                                ui *&temp_buffer) {
    VertexID u = order[depth];
    VertexID previous_bn = bn[depth][0];
    ui previous_index_id = idx_embedding[previous_bn];
    ui valid_candidates_count = 0;
    Edges& previous_edge = *edge_matrix[previous_bn][u];

    valid_candidates_count = previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui* previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

    memcpy(valid_candidate_index[depth], previous_candidates, valid_candidates_count * sizeof(ui));

    ui temp_count;
    for (ui i = 1; i < bn_cnt[depth]; ++i) {
        VertexID current_bn = bn[depth][i];
        Edges& current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count = current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui* current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, valid_candidate_index[depth], valid_candidates_count,
                        temp_buffer, temp_count);

        std::swap(temp_buffer, valid_candidate_index[depth]);
        valid_candidates_count = temp_count;
    }

    idx_count[depth] = valid_candidates_count;
}

bool
EvaluateQuery::GQLEngine(const Graph *data_graph, const Graph *query_graph, ui **candidates,ui *candidates_count,
                                   ui *order, uint64_t output_limit_num, uint64_t &call_cnt, 
                                   mpz_t embedding_cnt, int64_t& time_limit) {
    bool overtime = false;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    // Generate the bn.
    ui **bn;
    ui *bn_count;

    bn = new ui *[max_depth];
    for (ui i = 0; i < max_depth; ++i) {
        bn[i] = new ui[max_depth];
    }

    bn_count = new ui[max_depth];
    std::fill(bn_count, bn_count + max_depth, 0);

    std::vector<bool> visited_query_vertices(max_depth, false);
    visited_query_vertices[start_vertex] = true;
    for (ui i = 1; i < max_depth; ++i) {
        VertexID cur_vertex = order[i];
        ui nbr_cnt;
        const VertexID *nbrs = query_graph->getVertexNeighbors(cur_vertex, nbr_cnt);

        for (ui j = 0; j < nbr_cnt; ++j) {
            VertexID nbr = nbrs[j];

            if (visited_query_vertices[nbr]) {
                bn[i][bn_count[i]++] = nbr;
            }
        }

        visited_query_vertices[cur_vertex] = true;
    }

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    VertexID **valid_candidate;
    bool *visited_vertices;

    idx = new ui[max_depth];
    idx_count = new ui[max_depth];
    embedding = new ui[max_depth];
    visited_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
    valid_candidate = new ui *[max_depth];

    for (ui i = 0; i < max_depth; ++i) {
        VertexID cur_vertex = order[i];
        ui max_candidate_count = candidates_count[cur_vertex];
        valid_candidate[i] = new VertexID[max_candidate_count];
    }

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];
    std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
              valid_candidate[cur_depth]);

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
            embedding[u] = v;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            if (TimeOp::getClockNan() >= time_limit) {
                overtime = true;
                goto EXIT;
            }
            if ((ui)cur_depth == max_depth - 1) {
                mpz_add_ui(embedding_cnt, embedding_cnt, 1);
                visited_vertices[v] = false;
                if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
                    goto EXIT;
                }
            } else {
                call_cnt += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                        visited_vertices, bn, bn_count, order, candidates, candidates_count,
                                        query_graph);
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
    }

    // Release the buffer.
    EXIT:
    delete[] bn_count;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i) {
        delete[] bn[i];
        delete[] valid_candidate[i];
    }

    delete[] bn;
    delete[] valid_candidate;

    return overtime;
}

// 与CFL不同的是：这里没有pivot信息的使用，目前我(yujie)对pivot的理解是，VF3的parent信息
void EvaluateQuery::generateValidCandidates(const Graph *data_graph, ui depth, ui *embedding, ui *idx_count,
                                            ui **valid_candidate, bool *visited_vertices, ui **bn, ui *bn_cnt,
                                            ui *order, ui **candidates, ui *candidates_count,
                                            const Graph *query_graph) {
    VertexID u = order[depth];

    idx_count[depth] = 0;

    for (ui i = 0; i < candidates_count[u]; ++i) {
        VertexID v = candidates[u][i];

        if (!visited_vertices[v]) {
            bool valid = true;

            for (ui j = 0; j < bn_cnt[depth]; ++j) {
                VertexID u_bn = bn[depth][j];
                VertexID u_bn_v = embedding[u_bn];
#ifdef ELABELED_GRAPH
                LabelID elabel = query_graph->getEdgeLabel(u, u_bn, true);
                if (!data_graph->checkEdgeExistence(v, u_bn_v, elabel)) {
#else
                if (!data_graph->checkEdgeExistence(v, u_bn_v)) {
#endif
                    valid = false;
                    break;
                }
            }

            if (valid) {
                valid_candidate[depth][idx_count[depth]++] = v;
            }
        }
    }
}

bool
EvaluateQuery::QSIEngine(const Graph *data_graph, const Graph *query_graph, ui **candidates,ui *candidates_count,
                                   ui *order, ui *pivot, uint64_t output_limit_num, uint64_t &call_cnt,
                                   mpz_t embedding_cnt, int64_t& time_limit) {
    bool overtime = false;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    // Generate the bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, pivot, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    VertexID **valid_candidate;
    bool *visited_vertices;

    idx = new ui[max_depth];
    idx_count = new ui[max_depth];
    embedding = new ui[max_depth];
    visited_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
    valid_candidate = new ui *[max_depth];

    ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();
    for (ui i = 0; i < max_depth; ++i) {
        valid_candidate[i] = new VertexID[max_candidate_count];
    }

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];
    std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
              valid_candidate[cur_depth]);

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
            embedding[u] = v;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            if (TimeOp::getClockNan() >= time_limit) {
                overtime = true;
                goto EXIT;
            }
            if ((ui)cur_depth == max_depth - 1) {
                mpz_add_ui(embedding_cnt, embedding_cnt, 1);
                visited_vertices[v] = false;
                if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
                    goto EXIT;
                }
            } else {
                call_cnt += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(query_graph, data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                        visited_vertices, bn, bn_count, order, pivot);
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
    }

    // Release the buffer.
    EXIT:
    delete[] bn_count;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i) {
        delete[] bn[i];
        delete[] valid_candidate[i];
    }

    delete[] bn;
    delete[] valid_candidate;

    return overtime;
}

// 增加一个同label邻居数量的检查
void EvaluateQuery::generateValidCandidates(const Graph *query_graph, const Graph *data_graph, ui depth, ui *embedding,
                                            ui *idx_count, ui **valid_candidate, bool *visited_vertices, ui **bn,
                                            ui *bn_cnt,
                                            ui *order, ui *pivot) {
    VertexID u = order[depth];
    LabelID u_label = query_graph->getVertexLabel(u);
    ui u_degree = query_graph->getVertexDegree(u);

    idx_count[depth] = 0;

    VertexID p = embedding[pivot[depth]];
    ui nbr_cnt;
    auto nbrs = data_graph->getVertexNeighbors(p, nbr_cnt);
#ifdef ELABELED_GRAPH
    auto uqlabel = query_graph->getEdgeLabel(u, pivot[depth], true);
#endif

    for (ui i = 0; i < nbr_cnt; ++i) {
        VertexID v = nbrs[i];
#ifdef ELABELED_GRAPH
        if (data_graph->getEdgeLabel(p, nbrs[i], true) != uqlabel) continue;
#endif
        if (!visited_vertices[v] && u_label == data_graph->getVertexLabel(v) &&
            u_degree <= data_graph->getVertexDegree(v)) {
            bool valid = true;

            for (ui j = 0; j < bn_cnt[depth]; ++j) {
                VertexID u_bn = bn[depth][j];
                VertexID u_bn_v = embedding[u_bn];
#ifdef ELABELED_GRAPH
                LabelID elabel = query_graph->getEdgeLabel(u, u_bn, true);
                if (!data_graph->checkEdgeExistence(v, u_bn_v, elabel)) {
#else
                if (!data_graph->checkEdgeExistence(v, u_bn_v)) {
#endif
                    valid = false;
                    break;
                }
            }

            if (valid) {
                valid_candidate[depth][idx_count[depth]++] = v;
            }
        }
    }
}

bool
EvaluateQuery::DPisoEngine(const Graph *data_graph, const Graph *query_graph, TreeNode *tree,
                                        Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                        ui **weight_array, ui *order, uint64_t output_limit_num,
                                        uint64_t &call_cnt, mpz_t embedding_cnt, int64_t& time_limit) {
    ui max_depth = query_graph->getVerticesCount();

    ui *extendable = new ui[max_depth];
    for (ui i = 0; i < max_depth; ++i) {
        extendable[i] = tree[i].bn_count_;
    }

    // Generate backward neighbors.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    // Evaluate the query.
    bool overtime = false;
    int cur_depth = 0;

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, tree, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;  // v->u used for conflict
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    VertexID start_vertex = order[0];
    std::vector<dpiso_min_pq> vec_rank_queue;

    for (ui i = 0; i < candidates_count[start_vertex]; ++i) {
        VertexID v = candidates[start_vertex][i];
        embedding[start_vertex] = v;
        idx_embedding[start_vertex] = i;
        visited_vertices[v] = true;

#ifdef ENABLE_FAILING_SET
        reverse_embedding[v] = start_vertex;
#endif
        vec_rank_queue.emplace_back(dpiso_min_pq(extendable_vertex_compare));
        updateExtendableVertex(idx_embedding, idx_count, valid_candidate_idx, edge_matrix, temp_buffer, weight_array,
                               tree, start_vertex, extendable,
                               vec_rank_queue, query_graph);

        VertexID u = vec_rank_queue.back().top().first.first;
        vec_rank_queue.back().pop();

#ifdef ENABLE_FAILING_SET
        if (idx_count[u] == 0) {
            vec_failing_set[cur_depth] = ancestors[u];
        } else {
            vec_failing_set[cur_depth].reset();
        }
#endif

        call_cnt += 1;
        cur_depth += 1;
        order[cur_depth] = u;
        idx[u] = 0;
        while (cur_depth > 0) { // 回溯到第一层的时候结束
            while (idx[u] < idx_count[u]) { // 当前节点的candidate结束的时候，回溯
                ui valid_idx = valid_candidate_idx[u][idx[u]];
                v = candidates[u][valid_idx];

                if (visited_vertices[v]) {
                    idx[u] += 1;
#ifdef ENABLE_FAILING_SET // case 1:发生了冲突，FM = anc(u) ∪ anc(u′)
                    vec_failing_set[cur_depth] = ancestors[u];
                    vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                    // 上一层的点的failing set 与其子节点的failing set 取并集
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                    continue;
                }
                embedding[u] = v;
                idx_embedding[u] = valid_idx;
                visited_vertices[v] = true;
                idx[u] += 1;

#ifdef ENABLE_FAILING_SET
                reverse_embedding[v] = u;
#endif
                if (TimeOp::getClockNan() >= time_limit) {
                    overtime = true;
                    goto EXIT;
                }

                if ((ui)cur_depth == max_depth - 1) {
                    mpz_add_ui(embedding_cnt, embedding_cnt, 1);
                    visited_vertices[v] = false;
#ifdef ENABLE_FAILING_SET // case 3:匹配成功， 失败集置空（描述有问题，实际上应该是失败集中放入所有点）
                    reverse_embedding.erase(embedding[u]);
                    vec_failing_set[cur_depth].set();
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

#endif
                    if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
                        goto EXIT;
                    }
                } else { // 还有没有匹配的点，则向下走一层（匹配下一个点）
                    call_cnt += 1;
                    cur_depth += 1;
                    vec_rank_queue.emplace_back(vec_rank_queue.back());
                    updateExtendableVertex(idx_embedding, idx_count, valid_candidate_idx, edge_matrix, temp_buffer,
                                           weight_array, tree, u, extendable,
                                           vec_rank_queue, query_graph);

                    u = vec_rank_queue.back().top().first.first;
                    vec_rank_queue.back().pop();
                    idx[u] = 0;
                    order[cur_depth] = u;

#ifdef ENABLE_FAILING_SET
                    if (idx_count[u] == 0) { // 如果下一个点的候选集是空，那么直接判定失败
                        // case 2:因为候选集为空而失败  FM = anc(u)
                        vec_failing_set[cur_depth - 1] = ancestors[u];
                    } else { // 如果不是空，失败集信息清空
                        vec_failing_set[cur_depth - 1].reset();
                    }
#endif
                }
            }

            // 回溯
            cur_depth -= 1;
            vec_rank_queue.pop_back();
            u = order[cur_depth];
            visited_vertices[embedding[u]] = false;
            restoreExtendableVertex(tree, u, extendable);
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0) {
                if (!vec_failing_set[cur_depth].test(u)) {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    // 失败集那么多那么多的操作，实际剪枝就是这一行
                    // 如果下一个点出现没有出现在失败集中，那么下一个点不是失败原因，直接跳过
                    // 相对应这里的操作是，把candidate的索引指针 idx[u] ，指向末尾
                    idx[u] = idx_count[u];
                } else {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
        }
    }

    // Release the buffer.
    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    return overtime;
}

void EvaluateQuery::updateExtendableVertex(ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                           Edges ***edge_matrix, ui *&temp_buffer, ui **weight_array,
                                           TreeNode *tree, VertexID mapped_vertex, ui *extendable,
                                           std::vector<dpiso_min_pq> &vec_rank_queue, const Graph *query_graph) {
    TreeNode &node = tree[mapped_vertex];
    for (ui i = 0; i < node.fn_count_; ++i) {
        VertexID u = node.fn_[i];
        extendable[u] -= 1;
        if (extendable[u] == 0) {
            generateValidCandidateIndex(u, idx_embedding, idx_count, valid_candidate_index[u], edge_matrix, tree[u].bn_,
                                        tree[u].bn_count_, temp_buffer);

            ui weight = 0;
            for (ui j = 0; j < idx_count[u]; ++j) {
                ui idx = valid_candidate_index[u][j];
                weight += weight_array[u][idx];
            }
            vec_rank_queue.back().emplace(std::make_pair(std::make_pair(u, query_graph->getVertexDegree(u)), weight));
        }
    }
}

void EvaluateQuery::restoreExtendableVertex(TreeNode *tree, VertexID unmapped_vertex, ui *extendable) {
    TreeNode &node = tree[unmapped_vertex];
    for (ui i = 0; i < node.fn_count_; ++i) {
        VertexID u = node.fn_[i];
        extendable[u] += 1;
    }
}

void
EvaluateQuery::generateValidCandidateIndex(VertexID u, ui *idx_embedding, ui *idx_count, ui *&valid_candidate_index,
                                           Edges ***edge_matrix, ui *bn, ui bn_cnt, ui *&temp_buffer) {
    VertexID previous_bn = bn[0];
    Edges &previous_edge = *edge_matrix[previous_bn][u];
    ui previous_index_id = idx_embedding[previous_bn];

    ui previous_candidates_count =
            previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui *previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

    ui valid_candidates_count = 0;
    for (ui i = 0; i < previous_candidates_count; ++i) {
        valid_candidate_index[valid_candidates_count++] = previous_candidates[i];
    }

    ui temp_count;
    for (ui i = 1; i < bn_cnt; ++i) {
        VertexID current_bn = bn[i];
        Edges &current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count =
                current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, valid_candidate_index,
                                                  valid_candidates_count,
                                                  temp_buffer, temp_count);

        std::swap(temp_buffer, valid_candidate_index);
        valid_candidates_count = temp_count;
    }

    idx_count[u] = valid_candidates_count;
}

void EvaluateQuery::computeAncestor(const Graph *query_graph, TreeNode *tree, VertexID *order,
                                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors) {
    ui q_num = query_graph->getVerticesCount();
    ancestors.resize(q_num);

    // Compute the ancestor in the top-down order.
    for (ui i = 0; i < q_num; ++i) {
        VertexID u = order[i];
        TreeNode &u_node = tree[u];
        ancestors[u].set(u);
        for (ui j = 0; j < u_node.bn_count_; ++j) {
            VertexID u_bn = u_node.bn_[j];
            ancestors[u] |= ancestors[u_bn];
        }
    }
}

bool
EvaluateQuery::CECIEngine(const Graph *data_graph, const Graph *query_graph, TreeNode *tree, ui **candidates,
                                ui *candidates_count,
                                std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                                ui *order, uint64_t &output_limit_num, uint64_t &call_cnt,
                                mpz_t embedding_cnt, int64_t& time_limit) {

    ui max_depth = query_graph->getVerticesCount();
    ui data_vertices_count = data_graph->getVerticesCount();
    ui max_valid_candidates_count = 0;
    for (ui i = 0; i < max_depth; ++i) {
        if (candidates_count[i] > max_valid_candidates_count) {
            max_valid_candidates_count = candidates_count[i];
        }
    }
    // Allocate the memory buffer.
    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui *embedding = new ui[max_depth];
    ui *temp_buffer = new ui[max_valid_candidates_count];
    ui **valid_candidates = new ui *[max_depth];
    for (ui i = 0; i < max_depth; ++i) {
        valid_candidates[i] = new ui[max_valid_candidates_count];
    }
    bool *visited_vertices = new bool[data_vertices_count];
    std::fill(visited_vertices, visited_vertices + data_vertices_count, false);

    // Evaluate the query.
    bool overtime = false;
    int cur_depth = 0;
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidates[cur_depth][i] = candidates[start_vertex][i];
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidates[cur_depth][idx[cur_depth]];
            idx[cur_depth] += 1;

            if (visited_vertices[v]) {
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }

            embedding[u] = v;
            visited_vertices[v] = true;

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif
            if (TimeOp::getClockNan() >= time_limit) {
                overtime = true;
                goto EXIT;
            }
            if ((ui)cur_depth == max_depth - 1) {
                mpz_add_ui(embedding_cnt, embedding_cnt, 1);
                visited_vertices[v] = false;
#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
                    goto EXIT;
                }
            } else {
                call_cnt += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(cur_depth, embedding, idx_count, valid_candidates, order, temp_buffer, tree,
                                        TE_Candidates,
                                        NTE_Candidates);
#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0) {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                } else {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0) {
                if (!vec_failing_set[cur_depth].test(u)) {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                } else {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;
        }
    }

    // Release the buffer.
    EXIT:
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] temp_buffer;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i) {
        delete[] valid_candidates[i];
    }
    delete[] valid_candidates;

    return overtime;
}

void EvaluateQuery::generateValidCandidates(ui depth, ui *embedding, ui *idx_count, ui **valid_candidates, ui *order,
                                            ui *&temp_buffer, TreeNode *tree,
                                            std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates) {

    VertexID u = order[depth];
    idx_count[depth] = 0;
    ui valid_candidates_count = 0;
    {
        VertexID u_p = tree[u].parent_;
        VertexID v_p = embedding[u_p];

        auto iter = TE_Candidates[u].find(v_p);
        if (iter == TE_Candidates[u].end() || iter->second.empty()) {
            return;
        }

        valid_candidates_count = iter->second.size();
        VertexID *v_p_nbrs = iter->second.data();

        for (ui i = 0; i < valid_candidates_count; ++i) {
            valid_candidates[depth][i] = v_p_nbrs[i];
        }
    }
    ui temp_count;
    for (ui i = 0; i < tree[u].bn_count_; ++i) {
        VertexID u_p = tree[u].bn_[i];
        VertexID v_p = embedding[u_p];

        auto iter = NTE_Candidates[u][u_p].find(v_p);
        if (iter == NTE_Candidates[u][u_p].end() || iter->second.empty()) {
            return;
        }

        ui current_candidates_count = iter->second.size();
        ui *current_candidates = iter->second.data();

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count,
                                                  valid_candidates[depth], valid_candidates_count,
                                                  temp_buffer, temp_count);

        std::swap(temp_buffer, valid_candidates[depth]);
        valid_candidates_count = temp_count;
    }

    idx_count[depth] = valid_candidates_count;
}

void EvaluateQuery::computeAncestor(const Graph *query_graph, ui **bn, ui *bn_cnt, VertexID *order,
                                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors) {
    ui q_num = query_graph->getVerticesCount();
    ancestors.resize(q_num);

    // Compute the ancestor in the top-down order.
    for (ui i = 0; i < q_num; ++i) {
        VertexID u = order[i];
        ancestors[u].set(u);
        for (ui j = 0; j < bn_cnt[i]; ++j) {
            VertexID u_bn = bn[i][j];
            ancestors[u] |= ancestors[u_bn];
        }
    }
}

void EvaluateQuery::computeAncestor(const Graph *query_graph, VertexID *order,
                                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors) {
    ui q_num = query_graph->getVerticesCount();
    ancestors.resize(q_num);

    // Compute the ancestor in the top-down order.
    for (ui i = 0; i < q_num; ++i) {
        VertexID u = order[i];
        ancestors[u].set(u);
        for (ui j = 0; j < i; ++j) {
            VertexID u_bn = order[j];
            if (query_graph->checkEdgeExistence(u, u_bn)) {
                ancestors[u] |= ancestors[u_bn];
            }
        }
    }
}

bool
EvaluateQuery::VF3Engine(const Graph *data_graph, const Graph *query_graph, ui *order,
                               ui *pivot, uint64_t output_limit_num, uint64_t &call_cnt,
                               mpz_t embedding_cnt, int64_t& time_limit) {
    ui q_num = query_graph->getVerticesCount();
    ui query_labels_num = query_graph->getLabelsCount();
    ui query_max_label_fre = query_graph->getGraphMaxLabelFrequency();
    ui d_num = data_graph->getVerticesCount();
    ui data_labels_num = data_graph->getLabelsCount();
    ui data_max_label_fre = data_graph->getGraphMaxLabelFrequency();

    // 计算Feasibility sets: L, V
    // query_feaibility
    ui ***query_feasibility = NULL;
    ui **query_feasibility_count = NULL;
    query_feasibility = new ui**[q_num + 1]; // 这里取+1，参考原文
    query_feasibility_count = new ui*[q_num + 1];
    for (ui i = 0; i <= q_num; i++){
        query_feasibility[i] = new ui*[query_labels_num];
        query_feasibility_count[i] = new ui[query_labels_num];
        memset(query_feasibility_count[i], 0, query_labels_num*sizeof(ui));
        for (ui j = 0; j < query_labels_num; j++){
            query_feasibility[i][j] = new ui[query_max_label_fre];
        }
    }
    // 计算query_feasibility
    generateFeasibility(query_graph, order, query_feasibility, query_feasibility_count);

    // data_feasibility
    ui ***data_feasibility = NULL;
    ui **data_feasibility_count = NULL;
    data_feasibility = new ui**[q_num + 1];
    data_feasibility_count = new ui*[q_num + 1];
    for (ui i = 0; i <= q_num; i++){
        data_feasibility[i] = new ui*[data_labels_num];
        data_feasibility_count[i] = new ui[data_labels_num];
        memset(data_feasibility_count[i], 0, data_labels_num*sizeof(ui));
        for (ui j = 0; j < data_labels_num; j++){
            data_feasibility[i][j] = new ui[data_max_label_fre];
        }
    }

    // 匹配
    ui max_depth = query_graph->getVerticesCount();
    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui *embedding = new ui[max_depth];
    ui depth = 0;
    ui **valid_candidates = new ui *[max_depth];
    for (ui i = 0; i < max_depth; ++i) {
        valid_candidates[i] = new ui[data_max_label_fre];
    }
    bool *visited_vertices = new bool[d_num];
    memset(visited_vertices, 0, d_num*sizeof(bool));
    bool overtime = exploreVF3Backtrack(data_graph, query_graph, order, pivot, output_limit_num, call_cnt, embedding_cnt,
                        depth, max_depth, idx, idx_count, embedding, valid_candidates, visited_vertices,
                        query_feasibility, query_feasibility_count, data_feasibility, data_feasibility_count, time_limit);

    for (ui i = 0; i <= q_num; i++){
        for (ui j = 0; j < query_labels_num; j++){
            delete[] query_feasibility[i][j];
        }
        delete[] query_feasibility[i];
        delete[] query_feasibility_count[i];
    }
    delete[] query_feasibility;
    delete[] query_feasibility_count;
    for (ui i = 0; i <= q_num; i++){
        for (ui j = 0; j < data_labels_num; j++){
            delete[] data_feasibility[i][j];
        }
        delete[] data_feasibility[i];
        delete[] data_feasibility_count[i];
    }
    delete[] data_feasibility;
    delete[] data_feasibility_count;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    for (ui i = 0; i < max_depth; i++) {
        delete[] valid_candidates[i];
    }
    delete[] valid_candidates;
    delete[] visited_vertices;
    return overtime;
}

// 给定一个点，一张图，一组已匹配点，计算该点的feasibility set
void EvaluateQuery::generateFeasibility(const Graph *graph, ui v, bool* matched, ui level, ui*** feasibility, ui**feasibility_count) {
    // 初始化当前点的 feasibility set 信息：把上一个点的信息复制下来（刨去当前点）
    ui labels_num = graph->getLabelsCount();
    for (ui j = 0; j < labels_num; j++) {
        feasibility_count[level][j] = 0;
        for (ui k = 0; k < feasibility_count[level - 1][j]; k++) {
            if (feasibility[level - 1][j][k] != v) {
                feasibility[level][j][feasibility_count[level][j]++] = feasibility[level - 1][j][k];
            }
        }
    }
    ui nbrs_cnt;
    const ui *nbrs = graph->getVertexNeighbors(v, nbrs_cnt);
    for (ui j = 0; j < nbrs_cnt; j++) {
        if (!matched[nbrs[j]]) {
            // TODO, 先简单去个重试试，比较高级的去重需要更换数据结构
            ui label = graph->getVertexLabel(nbrs[j]);
            bool f_add = true;
            for (ui k  = 0; k < feasibility_count[level][label]; k++) {
                if (nbrs[j] == feasibility[level][label][k]) {
                    f_add = false;
                    break;
                }
            }
            if (f_add == true) {
                ui count = feasibility_count[level][label]++;
                feasibility[level][label][count] = nbrs[j];
            }
        }
    }
}

// 计算query_graph每个点的feasibility set
void EvaluateQuery::generateFeasibility(const Graph *query_graph, const ui *order, ui*** feasibility, ui**feasibility_count) {
    ui vertices_num = query_graph->getVerticesCount();

    bool* query_matched = new bool[vertices_num];
    memset(query_matched, 0, vertices_num*sizeof(bool));
    
    // 计算feasiblity信息，这里只需要记录连接的点即可，无向图
    // 最后一个点的信息不需要算，肯定是空的
    for (ui i = 1; i < vertices_num; i++) {
        ui u = order[i - 1];
        query_matched[u] = true;
        generateFeasibility(query_graph, u, query_matched, i, feasibility, feasibility_count);
    }
    delete[] query_matched;
}

bool EvaluateQuery::exploreVF3Backtrack(const Graph *data_graph, const Graph *query_graph,
                            ui *order, ui *pivot, uint64_t output_limit_num, uint64_t &call_cnt,
                            mpz_t embedding_cnt, ui depth, ui max_depth, ui*idx, ui*idx_count,
                            ui *embedding, ui **valid_candidates, bool* visited_vertices,
                            ui*** query_feasibility, ui**query_feasibility_count,
                            ui*** data_feasibility, ui**data_feasibility_count, int64_t& time_limit) {
    if (TimeOp::getClockNan() >= time_limit) {
        return  true;
    }
    if (depth == max_depth){
        // no checks so far
        mpz_add_ui(embedding_cnt, embedding_cnt, 1);
        if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
            for (ui i = 0; i < depth; i++){
                idx[i] = idx_count[i];
            }
        }
        return false;
    }
    ui u = order[depth];
    // 当前节点的 candidate set 的计算
    const ui*candidates;
    ui candidate_count;
    idx[depth] = 0;
    idx_count[depth] = 0;
    if (pivot[depth] == (ui)-1) {
        ui query_label = query_graph->getVertexLabel(u);
        candidates = data_graph->getVerticesByLabel(query_label, candidate_count);
        for (ui i = 0; i < candidate_count; i++) {
            if (!visited_vertices[candidates[i]]) {
                valid_candidates[depth][idx_count[depth]++] = candidates[i];
            }
        }
    } else {
        ui v = embedding[pivot[depth]];
        candidates = data_graph->getVertexNeighbors(v, candidate_count);
        ui query_label = query_graph->getVertexLabel(u);
        for (ui i = 0; i < candidate_count; i++) {
            ui data_label = data_graph->getVertexLabel(candidates[i]);
            if (!visited_vertices[candidates[i]] && data_label == query_label) {
                valid_candidates[depth][idx_count[depth]++] = candidates[i];
            }
        }
    }
    // 遍历当前节点的 valid_candidates
    while (idx[depth] < idx_count[depth]) {
        if (isFeasibility(data_graph, query_graph, depth, order[depth], valid_candidates[depth][idx[depth]],
                            embedding, order, query_feasibility, query_feasibility_count,
                            data_feasibility, data_feasibility_count) == true) {
            embedding[order[depth]] = valid_candidates[depth][idx[depth]];
            visited_vertices[embedding[order[depth]]] = true;
            call_cnt++;
            if (exploreVF3Backtrack(data_graph, query_graph, order, pivot, output_limit_num, call_cnt, embedding_cnt,
                        depth + 1, max_depth, idx, idx_count, embedding, valid_candidates, visited_vertices,
                        query_feasibility, query_feasibility_count, data_feasibility, data_feasibility_count, time_limit) == true) {
                return true;
            }
            visited_vertices[embedding[order[depth]]] = false;
        }
        idx[depth]++;
    }
    return false;
}

bool EvaluateQuery::isFeasibility(const Graph *data_graph, const Graph *query_graph, ui depth, ui cur_u, ui cur_v,
                    ui *embedding, ui* order, ui*** query_feasibility, ui**query_feasibility_count,
                    ui*** data_feasibility, ui**data_feasibility_count) {
    ui d_num = data_graph->getVerticesCount();
    ui q_num = query_graph->getVerticesCount();
    ui data_labels_num = data_graph->getLabelsCount();
    ui query_labels_num = query_graph->getLabelsCount();
    ui labels_num = query_labels_num;
    bool result = true;
    // 计算已匹配点 matched
    bool *data_matched = new bool[d_num];
    bool *query_matched = new bool[q_num];
    memset(data_matched, false, d_num*sizeof(bool));
    memset(query_matched, false, q_num*sizeof(bool));
    for (ui i = 0; i < depth; i++) {
        data_matched[embedding[order[i]]] = true;
        query_matched[order[i]] = true;
    }
    // 求数据图的feasibility set
    generateFeasibility(data_graph, cur_v, data_matched, depth + 1, data_feasibility, data_feasibility_count);

    // 计算feasibility set matched
    bool *data_fmatched = new bool[d_num];
    bool *query_fmatched = new bool[q_num];
    memset(data_fmatched, false, d_num*sizeof(bool));
    memset(query_fmatched, false, q_num*sizeof(bool));
    for (ui i = 0; i < query_labels_num; i++) {
        for (ui j = 0; j < query_feasibility_count[depth+1][i]; j++) {
            query_fmatched[query_feasibility[depth+1][i][j]] = true;
        }
    }
    for (ui i = 0; i < data_labels_num; i++) {
        for (ui j = 0; j < data_feasibility_count[depth+1][i]; j++) {
            data_fmatched[data_feasibility[depth+1][i][j]] = true;
        }
    }

    // feasibility rules 信息处理（处理两点的neighbors）
    ui unbrs_count = 0 ;
    auto unbrs = query_graph->getVertexNeighbors(cur_u, unbrs_count);
#ifdef ELABELED_GRAPH
    auto elabels = query_graph->getVertexEdgeLabels(cur_u, unbrs_count);
#endif
    ui vnbrs_count = 0;
    const ui *vnbrs = data_graph->getVertexNeighbors(cur_v, vnbrs_count);
    if (vnbrs_count < unbrs_count) {
        result = false;
    } else {
        // 用于计算不同label邻居的数量(in feasibility)
        ui *unbrs_flabeled_num = new ui[query_labels_num];
        ui *vnbrs_flabeled_num = new ui[data_labels_num];
        memset(unbrs_flabeled_num, 0, query_labels_num*sizeof(ui));
        memset(vnbrs_flabeled_num, 0, data_labels_num*sizeof(ui));
        // 用于计算不同label邻居的数量(not in feasibility & not matched)
        ui *unbrs_labeled_num = new ui[query_labels_num];
        ui *vnbrs_labeled_num = new ui[data_labels_num];
        memset(unbrs_labeled_num, 0, query_labels_num*sizeof(ui));
        memset(vnbrs_labeled_num, 0, data_labels_num*sizeof(ui));

        for (ui i = 0; i < unbrs_count; i++) {
            ui unbr = unbrs[i];
            ui unbr_label = query_graph->getVertexLabel(unbr);
            
            if (query_matched[unbr] == true) {
#ifdef ELABELED_GRAPH
                if (!data_graph->checkEdgeExistence(cur_v, embedding[unbr], elabels[i]))
#else
                if (!data_graph->checkEdgeExistence(cur_v, embedding[unbr]))
#endif
                    result = false;
            } else if (query_fmatched[unbr] == true) {
                unbrs_flabeled_num[unbr_label]++;
            } else {
                unbrs_labeled_num[unbr_label]++;
            }
        }

        for (ui i = 0; i < vnbrs_count; i++) {
            ui vnbr = vnbrs[i];
            ui vnbr_label = data_graph->getVertexLabel(vnbr);
            
            if (data_matched[vnbr] == true) {
                ; // do nothing;
            } else if (data_fmatched[vnbr] == true) {
                vnbrs_flabeled_num[vnbr_label]++;
            } else {
                vnbrs_labeled_num[vnbr_label]++;
            }
        }

        for (ui i = 0; i < labels_num; i++) {
            // for debug and believe the compiler
            if (vnbrs_labeled_num[i] < unbrs_labeled_num[i]) {
                result = false;
                break;
            }
            if (vnbrs_flabeled_num[i] < unbrs_flabeled_num[i]){
                result = false;
                break;
            }
        }
        delete[] unbrs_flabeled_num;
        delete[] vnbrs_flabeled_num;
        delete[] unbrs_labeled_num;
        delete[] vnbrs_labeled_num;
    }
    delete[] data_matched;
    delete[] query_matched;
    delete[] data_fmatched;
    delete[] query_fmatched;

    return result;
}

bool
EvaluateQuery::VEQEngine(const Graph *data_graph, const Graph *query_graph, TreeNode *tree,
                                    Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                    uint64_t output_limit_num, uint64_t &call_cnt,
                                    mpz_t embedding_cnt, int64_t& time_limit) {
    // 首先计算NEC，用于动态生成匹配顺序
    ui q_num = query_graph->getVerticesCount();
    ui d_num = data_graph->getVerticesCount();
    ui** nec = new ui*[q_num];
    memset(nec, 0, sizeof(ui*)*q_num);
    computeNEC(query_graph, nec);

    // 初始化辅助结构
    ui depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    ui *embedding = new ui[q_num]; // u->v
    ui* order = new ui[max_depth];      // depth->u 当前一层的u是谁
    bool overtime = false;               // 是否超时
    ui *extendable = new ui[q_num]; // 计算query_graph上的可扩展点
    for (ui i = 0; i < q_num; ++i) {
        extendable[i] = tree[i].bn_count_;
    }
    ui *idx = new ui[max_depth];                   // depth->cans_idx
    ui *idx_count = new ui[max_depth];             // depth->cans_count
    ui** valid_cans = new ui*[q_num]; // u->valid_cans
    for (ui i = 0; i < q_num; i++) {
        valid_cans[i] = new ui[candidates_count[i]];
    }
    ui* valid_cans_count = new ui[q_num]; // u->valid_cans_count
    bool* visited_u = new bool[q_num]; // vertexID->bool(某个点是否访问过)(query)
    memset(visited_u, 0, q_num*sizeof(bool));
    bool *visited_vertices = new bool[d_num]; // VertexID->bool(某个点是否访问过)(data)
    memset(visited_vertices, 0, d_num*sizeof(bool));

#ifdef ENABLE_EQUIVALENT_SET
    ui** TM = new ui*[q_num];   // 用于计算每个节点valid_cans中每个节点对应的子树中生成的匹配数
    for (ui i = 0; i < q_num; i++) {
        TM[i] = new ui[candidates_count[i] + 1]; // 最后一个位置用来放当前层的总数，每个点valid_cans最多有can_cnt个
    }
    memset(TM[0], 0, sizeof(ui)*(candidates_count[0] + 1));
     // vec_index和vec_set用于计算共享邻居的点
    std::vector<std::vector<ui>> vec_index(q_num);
    for (ui i = 0; i < q_num; i++) {
        vec_index[i].resize(candidates_count[i]);
        std::fill(vec_index[i].begin(), vec_index[i].end(), (ui)-1);
    }
    std::vector<std::vector<ui>> vec_set;  //预先分配的空间绝对不够
    // 这里vec_set是静态计算的，但是论文里写的是动态计算（之后再做考量吧）（TODO）
    overtime = computeNEC(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set, time_limit);
    if (overtime) return overtime;
    // pi_m_index和pi_m用于计算运行时的等价点
    std::vector<std::vector<ui>> pi_m_index(q_num);
    for (ui i = 0; i < q_num; i++) {
        pi_m_index[i].resize(candidates_count[i], (ui)-1);
    }
    std::vector<std::vector<ui>> pi_m;  // π_m
    ui* pi_m_count = new ui[max_depth]; //记录每一层的pi_m_index用到哪个数了
    pi_m_count[0] = 0;
    // dm_index和dm用于计算成功集信息
    std::vector<std::vector<ui>> dm(q_num);  // δ_m
    std::unordered_map<VertexID, VertexID> reverse_embedding;  // v->u used for conflict
    reverse_embedding.reserve(q_num); // 与embedding一样，最多有query_v_n个点
#endif

    // 第一个匹配点
    ui start_vertex = (ui)-1;
    for (ui i = 0; i < q_num; i++) {
        if (extendable[i] == 0) {
            start_vertex = i;
            break;
        }
    }
    assert(start_vertex != (ui)-1);
    order[depth] = start_vertex;
    visited_u[start_vertex] = true;
    // 计算第一个点的valid_cans和idx
    idx[depth] = 0;
    valid_cans_count[start_vertex] = candidates_count[start_vertex];
    idx_count[depth] = valid_cans_count[start_vertex];
    for (ui i = 0; i < candidates_count[start_vertex]; i++) {
        valid_cans[start_vertex][i] = candidates[start_vertex][i];
    }
#ifdef ENABLE_EQUIVALENT_SET
    std::fill(pi_m_index[start_vertex].begin(), pi_m_index[start_vertex].end(), (ui)-1);
#endif
    // 开始匹配
    while (true) {
        while (idx[depth] < idx_count[depth]) {
            if (TimeOp::getClockNan() >= time_limit) {
                overtime = true;
                goto EXIT;
            }
            // compute next u
            VertexID u = order[depth];
            VertexID v = valid_cans[u][idx[depth]];

            if (visited_vertices[v]) {
#ifdef ENABLE_EQUIVALENT_SET // 发生冲突， 处理产生冲突的节点
                TM[u][idx[depth]] = 0;
                VertexID con_u = reverse_embedding[v];
                ui con_v_index = 0, v_index = 0;
                for (; con_v_index < valid_cans_count[con_u]; con_v_index++) {
                    if (valid_cans[con_u][con_v_index] == v) break;
                }
                for (; v_index < candidates_count[u]; v_index++) {
                    if (candidates[u][v_index] == v) break;
                }
                assert(con_v_index < valid_cans_count[con_u]);
                assert(v_index < candidates_count[u]);
                // 求交
                auto& con_uv_idx = pi_m_index[con_u][con_v_index];
                for (ui i = 0; i < pi_m[con_uv_idx].size(); i++) {
                    bool f_in = false;
                    for (ui j = 0; j < vec_set[vec_index[u][v_index]].size(); j++) {
                        if (vec_set[vec_index[u][v_index]][j] == pi_m[con_uv_idx][i]) {
                            f_in = true;
                            break;
                        }
                    }
                    if (f_in == false) {
                        pi_m[con_uv_idx][i] = pi_m[con_uv_idx][pi_m[con_uv_idx].size() - 1];
                        pi_m[con_uv_idx].pop_back();
                        i--;
                    }
                }
#endif
                idx[depth]++;
                continue;
            }
#ifdef ENABLE_EQUIVALENT_SET // 如果当前节点是等价节点，跳过
            if (pi_m_index[u][idx[depth]] != (ui)-1) {
                // 一个小操作，如果是成功集的扩展，那么把导致成功的点，放到pi_m的第一位，需要记录匹配数
                VertexID equ_v = pi_m[pi_m_index[u][idx[depth]]][0];
                ui equ_v_index = 0;
                for (; equ_v_index < idx_count[depth]; equ_v_index++) {
                    if (equ_v == valid_cans[u][equ_v_index]) break;
                }
                assert(equ_v_index < idx_count[depth]);
                mpz_add_ui(embedding_cnt, embedding_cnt, TM[u][equ_v_index]);
                TM[u][candidates_count[u]] += TM[u][equ_v_index];
                idx[depth]++;
                continue;
            }
            reverse_embedding[v] = u;
#endif
            embedding[u] = v;
            visited_vertices[v] = true;
            ui cur_idx = idx[depth]++;
#ifdef ENABLE_EQUIVALENT_SET // 初始化 pi_m(u,v) 信息
            pi_m_index[u][cur_idx] = pi_m_count[depth]++;
            ui v_index = 0;
            for (; v_index < candidates_count[u]; v_index++) {
                if (candidates[u][v_index] == v) break;
            }
            auto& pi = vec_set[vec_index[u][v_index]];
            pi_m.push_back(pi); // 𝜋−𝑀(𝑢, 𝑣) ← 𝜋 (𝑢, 𝑣);
            assert(pi_m.size() == pi_m_count[depth]);
            dm[u].clear();                                  // 𝛿𝑀 (𝑢, 𝑣) ← ∅;
            for (ui i = 0; i < depth; i++) {                // for each ancestor
                if (TimeOp::getClockNan() >= time_limit) {
                    overtime = true;
                    goto EXIT;
                }
                bool va_in_pi = false, sec_empty = true;
                ui va_index = 0;
                VertexID ua = order[i];
                VertexID va = embedding[ua];
                for (; va_index < candidates_count[ua]; va_index++) {
                    if (candidates[ua][va_index] == va) break;
                }
                assert(va_index < candidates_count[ua]);
                auto& pi_a = vec_set[vec_index[ua][va_index]];
                for (ui j = 0; j < pi.size(); j++) { // 𝜋 (𝑢, 𝑣)
                    if (pi[j] == embedding[order[i]]) { // 𝑣𝑎 ∉ 𝜋 (𝑢, 𝑣)
                        va_in_pi = true;
                        break;
                    }
                    for (ui k = 0; k < pi_a.size(); k++) { // 𝜋 (𝑢𝑎,𝑣𝑎)∩𝜋 (𝑢, 𝑣) ≠ ∅
                        if (pi_a[k] == pi[j]) {
                            sec_empty = false;
                            break;
                        }
                    }
                    if (sec_empty == false) break;
                }
                if (va_in_pi == false && sec_empty == false) { // 𝛿𝑀 (𝑢𝑎,𝑣𝑎) ← 𝛿𝑀 (𝑢𝑎, 𝑣𝑎) ∪ 𝜋 (𝑢, 𝑣);
                    ui dm_size = dm[ua].size();
                    for (ui j = 0; j < pi.size(); j++) {
                        bool f_add = true;
                        for (ui k = 0; k < dm_size; k++) {
                            if (dm[ua][k] == pi[j]) {
                                f_add = false;
                                break;
                            }
                        }
                        if (f_add == true) {
                            dm[ua].push_back(pi[j]);
                        }
                    }
                }
            }
#endif
            if (depth == max_depth - 1) { // 匹配成功
                mpz_add_ui(embedding_cnt, embedding_cnt, 1);
                visited_vertices[v] = false;
                if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
                    goto EXIT;
                }
#ifdef ENABLE_EQUIVALENT_SET // 在最后一个点上匹配成功，将成功结果扩散到整层
                reverse_embedding.erase(v);
                TM[u][cur_idx] = 1;
                TM[u][candidates_count[u]]++;
                auto& uv_idx = pi_m_index[u][cur_idx];
                for (ui i = 0; i < pi_m[uv_idx].size(); i++) { // 𝜋𝑀 (𝑢, 𝑣) ← 𝜋−𝑀(𝑢, 𝑣) − 𝛿𝑀 (𝑢, 𝑣);
                    bool f_del = false;
                    for (ui j = 0; j < dm[u].size(); j++) {
                        if (pi_m[uv_idx][i] == dm[u][j]) {
                            f_del = true;
                            break;
                        }
                    }
                    if (f_del == true) {
                        pi_m[uv_idx][i] = pi_m[uv_idx][pi_m[uv_idx].size() - 1];
                        pi_m[uv_idx].pop_back();
                        i--;
                    }
                }
                for (ui i = 1; i < pi_m[uv_idx].size(); i++) { // foreach 𝑣′ ∈ 𝜋𝑀 (𝑢, 𝑣)
                    ui v_equ_index = 0;
                    for (; v_equ_index < valid_cans_count[u]; v_equ_index++) {
                        if (valid_cans[u][v_equ_index] == pi_m[uv_idx][i])
                            break;
                    }
                    assert(v_equ_index < valid_cans_count[u]);
                    pi_m_index[u][v_equ_index] = uv_idx;
                }
                // 将计算出的pi中导致成功的点，放到第一位上
                if (pi_m[uv_idx].size() > 1) {
                    auto pi_v_idx = std::find(pi_m[uv_idx].begin(), pi_m[uv_idx].end(), v);
                    assert(pi_v_idx != pi_m[uv_idx].end());
                    ui tmp = *pi_v_idx;
                    *pi_v_idx = pi_m[uv_idx][0];
                    pi_m[uv_idx][0] = tmp;
                }
#endif
            } else { // 匹配下一层
                call_cnt += 1;
                depth += 1;
                order[depth] = generateNextU(data_graph, query_graph, candidates, candidates_count, valid_cans,
                                             valid_cans_count, extendable, nec, depth, embedding,
                                             edge_matrix, visited_vertices, visited_u, order, tree);
                
                if (order[depth] == (ui)-1) {
                    break; // 对呀，这里用break就可以达到直接回溯的效果呢
                } else {
                    visited_u[order[depth]] = true;
                    idx[depth] = 0;
                    idx_count[depth] = valid_cans_count[order[depth]];
                }
#ifdef ENABLE_EQUIVALENT_SET // 初始化下一层节点的 等价集 信息
                memset(TM[order[depth]], 0, sizeof(ui)*(candidates_count[order[depth]] + 1));
                std::fill(pi_m_index[order[depth]].begin(), pi_m_index[order[depth]].end(), (ui)-1);
                pi_m_count[depth] = pi_m_count[depth - 1];
#endif
            }
        }
        // 回溯部分
        depth -= 1;
        if (depth == (ui)-1)
            break;
        VertexID u = order[depth];
        ui cur_idx = idx[depth] - 1;
        visited_vertices[embedding[u]] = false;
        restoreExtendableVertex(tree, u, extendable);
        if (order[depth + 1] != (ui)-1) { // 如果当前这一层没有找到合适的匹配点，那么会直接回溯
            VertexID last_u = order[depth + 1];
            visited_u[last_u] = false;
            if (nec[last_u] != NULL) (*(nec[last_u]))++;
        }
#ifdef ENABLE_EQUIVALENT_SET // 处理结束上一个节点，将结果扩散
        if (order[depth + 1] != (ui)-1) {
            TM[u][cur_idx] = TM[order[depth + 1]][candidates_count[order[depth + 1]]];
        } else {
            TM[u][cur_idx] = 0;
        }
        TM[u][candidates_count[u]] += TM[u][cur_idx];
        auto& uv_idx = pi_m_index[u][cur_idx];
        if (TM[u][cur_idx] != 0) {
            for (ui i = 0; i < pi_m[uv_idx].size(); i++) {
                bool f_del = false;
                for (ui j = 0; j < dm[u].size(); j++) {
                    if (pi_m[uv_idx][i] == dm[u][j]) {
                        f_del = true;
                        break;
                    }
                }
                if (f_del == true) {
                    pi_m[uv_idx][i] = pi_m[uv_idx][pi_m[uv_idx].size() - 1];
                    pi_m[uv_idx].pop_back();
                    i--;
                }
            }
            // 将计算出的pi中导致成功的点，放到第一位上
            if (pi_m[uv_idx].size() > 1) {
                auto pi_v_idx = std::find(pi_m[uv_idx].begin(), pi_m[uv_idx].end(), embedding[u]);
                assert(pi_v_idx != pi_m[uv_idx].end());
                ui tmp = *pi_v_idx;
                *pi_v_idx = pi_m[uv_idx][0];
                pi_m[uv_idx][0] = tmp;
            }
        }
        for (ui i = 1; i < pi_m[uv_idx].size(); i++) {
            ui v_equ_index = 0;
            for (; v_equ_index < valid_cans_count[u]; v_equ_index++) {
                if (valid_cans[u][v_equ_index] == pi_m[uv_idx][i])
                    break;
            }
            assert(v_equ_index < valid_cans_count[u]);
            pi_m_index[u][v_equ_index] = uv_idx;
        }
        pi_m.resize(pi_m_count[depth]);
#endif
    }

    EXIT:
    // 清理空间(部分空间先留着不请了，想到办法了再说)
    for (ui i = 0; i < q_num; i++) {
        // if (nec[i] != NULL) delete nec[i];
    }
    delete []nec;
    delete []embedding;
    delete []visited_u;
    delete []visited_vertices;
    delete []order;
    delete []extendable;
    delete []idx;
    delete []idx_count;
    for (ui i = 0; i < q_num; i++) {
        delete[] valid_cans[i];
    }
    delete []valid_cans;
    delete []valid_cans_count;
#ifdef ENABLE_EQUIVALENT_SET
for (ui i = 0; i < q_num; i++) {
        delete[] TM[i];
    }
    delete[] TM;
    delete[] pi_m_count;
#endif
    return overtime;
}

void EvaluateQuery::RestoreValidCans(const Graph *query_graph, const Graph *data_graph, bool* visited_u,
                                     VertexID last_u, VertexID last_v,
                                     std::vector<std::unordered_map<VertexID, ui>>& valid_cans) {
    ui last_unbrs_count;
    const ui* last_unbrs = query_graph->getVertexNeighbors(last_u, last_unbrs_count);
    ui last_vnbrs_count;
    const ui* last_vnbrs = data_graph->getVertexNeighbors(last_v, last_vnbrs_count);
    for (ui i = 0; i < last_unbrs_count; i++) {
        ui last_unbr = last_unbrs[i];
        if (visited_u[last_unbr] == true) continue;
        for (ui j = 0; j < last_vnbrs_count; j++) {
            auto vertex = valid_cans[last_unbr].find(last_vnbrs[j]);
            if (vertex != valid_cans[last_unbr].end()) {
                if (vertex->second == 1) {
                    valid_cans[last_unbr].erase(vertex);
                } else {
                    vertex->second--;
                }
            }
        }
    }
}

ui EvaluateQuery::generateNextU(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count,
                                ui**valid_cans, ui*valid_cans_count, ui* extendable,  ui** nec,
                                ui depth, ui* embedding, Edges ***edge_matrix, bool *visited_vertices,
                                bool *visited_u, ui *order, TreeNode* tree) {
    // 首先更新所有被影响到的点extenable点的valid_cans信息，用于计算NextU
    ui q_num = query_graph->getVerticesCount();
    ui cur_vertex = -1;
    TreeNode &node = tree[order[depth - 1]];
    for (ui i = 0; i < node.fn_count_; ++i) {
        VertexID u = node.fn_[i];
        extendable[u] -= 1;
        if (extendable[u] == 0) { // 当点前点的前向点(也就是指向其的点)结束后，这些点可以扩展
            ComputeValidCans(data_graph, query_graph, candidates, candidates_count, valid_cans,
                             valid_cans_count, embedding, u, visited_u);
        }
    }
    // 现在为止，valid_cans_count表示的就是|Cm(u)|，而valid_cans存的就是Cm(u)
    bool f_only_1 = true;
    // 先找满足|NEC(u)| >= |Um(u)|的点，同时记录degree == 1的点的情况
    for (ui i = 0; i < q_num; i++) {
        if (visited_u[i] == true || extendable[i] != 0) continue;
        ui u_count = 0;
        for (ui j = 0; j < valid_cans_count[i]; j++) {
            if (visited_vertices[valid_cans[i][j]] == false) u_count++;
        }
        if (nec[i] != NULL && *(nec[i]) >= u_count) {
            if (*(nec[i]) > u_count) {
                return (ui)-1;
            } else {
                (*(nec[i]))--;
                return i;
            }
        }
        if (nec[i] != NULL) {
            if(cur_vertex == (ui)-1) cur_vertex = i;
        } else if (f_only_1 == true) {
            f_only_1 = false;
        }
    }

    if (f_only_1 == false) { // 最后找|Um(u)|最小的degree != 1的点
        cur_vertex = (ui)-1;
        ui min_Um_u = (ui)-1;
        for (ui i = 0; i < q_num; i++) {
            if (visited_u[i] == true || nec[i] != 0 || extendable[i] != 0) continue;
            if (valid_cans_count[i] != 0 && min_Um_u > valid_cans_count[i]) {
                cur_vertex = i;
                min_Um_u = valid_cans_count[i];
            }
        }
    }

    if (cur_vertex != (ui)-1 && nec[cur_vertex] != NULL) {
        (*nec[cur_vertex])--;
    }

    return cur_vertex;
}

void EvaluateQuery::ComputeValidCans(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count,
                      ui**valid_cans, ui*valid_cans_count, ui* embedding, VertexID u, bool* visited_u) {
    ui unbrs_count;
    auto unbrs = query_graph->getVertexNeighbors(u, unbrs_count);
#ifdef ELABELED_GRAPH
    auto elabels = query_graph->getVertexEdgeLabels(u, unbrs_count);
#endif
    valid_cans_count[u] = 0;
    for (ui i = 0; i < candidates_count[u]; i++) {
        VertexID v = candidates[u][i];
        bool flag = true;
        for (ui j = 0; j < unbrs_count; j++) {
            if (visited_u[unbrs[j]] == true
#ifdef ELABELED_GRAPH
                && !data_graph->checkEdgeExistence(v, embedding[unbrs[j]], elabels[j])) {
#else
                && !data_graph->checkEdgeExistence(v, embedding[unbrs[j]])) {
#endif
                flag = false;
                break;
            }
        }
        if (flag == true) {
            valid_cans[u][valid_cans_count[u]++] = v;
        }
    }
}

// 计算candidates上的等价邻居
bool
EvaluateQuery::computeNEC(const Graph *query_graph, Edges ***edge_matrix, ui *candidates_count,
                               ui**candidates, std::vector<std::vector<ui>>& vec_index,
                               std::vector<std::vector<ui>>& vec_set, int64_t& time_limit) {
    std::vector<ui> tmp_vec;
    ui vec_count = 0;
    for (ui i = 0; i < vec_index.size(); i++) {
        tmp_vec.reserve(candidates_count[i]);
        ui unbrs_count;
        const ui *unbrs = query_graph->getVertexNeighbors(i, unbrs_count);
        for (ui j = 0; j < candidates_count[i]; j++) {
            if (vec_index[i][j] != (ui)-1)
                continue;
            vec_index[i][j] = vec_count++;
            tmp_vec.push_back(candidates[i][j]);
            for (ui k = j + 1; k < candidates_count[i]; k++) {
                if (TimeOp::getClockNan() >= time_limit) {
                    return true;
                }
                if (vec_index[i][k] != (ui)-1)
                    continue;
                // 这里判断邻居是否相同，并存入tmp_vec中
                bool equ = true;
                for (ui u1 = 0; u1 < unbrs_count; u1++) {
                    ui unbr = unbrs[u1];
                    const Edges* edges = edge_matrix[i][unbr];
                    if (edges->offset_[j+1]-edges->offset_[j] == 0) {
                        // 如果当前点的邻居的边数为 0，那么直接跳过
                        // 不应该放在任何一个等价集中
                        goto SKIP;
                    }
                    if (edges->offset_[j+1]-edges->offset_[j] != edges->offset_[k+1] - edges->offset_[k]) {
                        equ = false;
                        break;
                    }
                    // 这里默认边的顺序是按照candidates的顺序来的吧
                    // 减少一层for循环吧，实在是太深了
                    for (ui u2 = 0; u2 < edges->offset_[j+1] - edges->offset_[j]; u2++) {
                        if (edges->edge_[u2+edges->offset_[j]] != edges->edge_[u2+edges->offset_[k]]) {
                            equ = false;
                            break;
                        }
                    }
                    if (equ == false) break;
                }
                if (equ == true) {
                    tmp_vec.push_back(candidates[i][k]);
                    vec_index[i][k] = vec_index[i][j];
                }
            }
            SKIP:
            // 处理当前一批的vec
            vec_set.push_back(tmp_vec);
            tmp_vec.clear();
        }
    }
    return false;
}

// 用于计算query_graph上的度为 1 的等价点
void EvaluateQuery::computeNEC(const Graph *query_graph, ui** nec) {
    ui q_num = query_graph->getVerticesCount();
    ui* nec_tmp = new ui[q_num];
    ui flag = true; // 定义一个通用标识
    for (ui i = 0; i < q_num; i++) {
        if (query_graph->getVertexDegree(i) != 1 || nec[i] != NULL) {
            continue;
        }
        ui unbrs_num;
        auto unbrs = query_graph->getVertexNeighbors(i, unbrs_num);
#ifdef ELABELED_GRAPH
        auto u1elabels = query_graph->getVertexEdgeLabels(i, unbrs_num);
#endif
        ui u_label = query_graph->getVertexLabel(i);
        ui* nec_count = new ui(0);
        nec_tmp[(*nec_count)++] = i;
        for (ui j = i + 1; j < q_num; j++) {
            if (query_graph->getVertexDegree(i) != 1 || nec[i] != NULL) {
                continue;
            }
            // 比较label和nbr信息
            if (u_label != query_graph->getVertexLabel(j)) {
                continue;
            }
            ui u2_nbrs_num;
            auto u2_nbrs = query_graph->getVertexNeighbors(j, u2_nbrs_num);
#ifdef ELABELED_GRAPH
            auto u2elabels = query_graph->getVertexEdgeLabels(j, u2_nbrs_num);
#endif
            if (u2_nbrs_num != unbrs_num) {
                continue;
            }
            ui k1 = 0;
            for (k1 = 0; k1 < unbrs_num; k1++) {
                flag = false;
                for (ui k2 = 0; k2 < u2_nbrs_num; k2++) {
                    if (unbrs[k1] == u2_nbrs[k2]) {
#ifdef ELABELED_GRAPH
                        if (u1elabels[k1] == u2elabels[k2])
#endif
                            flag = true;
                        break;
                    }
                }
                if (flag == false) break; // 任意一个邻居不同，提前退出
            }
            if (k1 != unbrs_num) continue;
            nec_tmp[(*nec_count)++] = j;
        }
        // 记录找到的所有与i有关的nec
        for (ui j = 0; j < *nec_count; j++) {
            nec[nec_tmp[j]] = nec_count;
        }
    }
    delete []nec_tmp;
}

bool
EvaluateQuery::RMEngine(const Graph *query_graph, const Graph *data_graph, catalog*&storage, Edges ***edge_matrix,
                              ui **candidates, ui *candidates_count, ui *order, uint64_t output_limit_num, uint64_t &call_cnt,
                              mpz_t embedding_cnt, int64_t& time_limit) {
    // first construct the catalog info
#ifdef ELABELED_GRAPH
    if (storage != NULL) {
        delete storage;
        storage = NULL;
    }
#endif
    if (storage == NULL) {
        storage = new catalog(query_graph, data_graph);
        convertCans2Catalog(query_graph, candidates, edge_matrix, storage);
        storage->query_graph_->get2CoreSize();
        storage->data_graph_->getVerticesCount();
    }
    convert_to_encoded_relation(storage, order);
#ifdef SPARSE_BITMAP
    convert_encoded_relation_to_sparse_bitmap(storage, order);
#endif
    std::vector<ui> vorder;
    vorder.reserve(query_graph->getVerticesCount());
    vorder.insert(vorder.end(), order, order + query_graph->getVerticesCount());
    auto tree = execution_tree_generator::generate_single_node_execution_tree(vorder);
    size_t result_cnt = 0;
    auto overtime = tree->execute(*storage, output_limit_num, call_cnt, result_cnt, time_limit);
    mpz_init_set_ui(embedding_cnt, result_cnt);
    return overtime;
}

void EvaluateQuery::convertCans2Catalog(const Graph *query_graph, ui **candidates, Edges ***edge_matrix, catalog *storage) {
    // fill edge_realation_
    for (ui u = 0; u < storage->num_sets_; u++) {
        ui unbrs_cnt;
        const ui* unbrs = query_graph->getVertexNeighbors(u, unbrs_cnt);
        for (ui i = 0; i < unbrs_cnt; i++) {
            ui v = unbrs[i];
            // only add edges (src < dst)
            if (u > v) continue;
            std::vector<edge> tmp_edges;
            auto& edges = *edge_matrix[u][v];
            auto& relation = storage->edge_relations_[u][v];
            for (ui j = 0; j < edges.vertex_count_; j++) {
                ui src = candidates[u][j];
                for (ui k = edges.offset_[j]; k < edges.offset_[j+1]; k++) {
                    ui dst = candidates[v][edges.edge_[k]];
                    tmp_edges.push_back(std::move(edge(src, dst)));
                }
            }
            relation.size_ = tmp_edges.size();
            relation.edges_ = new edge[relation.size_];
            memcpy(relation.edges_, tmp_edges.data(), sizeof(edge) * relation.size_);
        }
    }
}

void EvaluateQuery::convert_to_encoded_relation(catalog *storage, ui*order) {
    auto& query_graph = storage->query_graph_;
    uint32_t core_vertices_cnt = query_graph->get2CoreSize();
    auto max_vertex_id = storage->data_graph_->getVerticesCount();

    auto projection_operator = new projection(max_vertex_id);
    for (uint32_t i = 0; i < core_vertices_cnt || i == 0; ++i) {
        uint32_t u = order[i];
        uint32_t nbr_cnt;
        const uint32_t* nbrs = query_graph->getVertexNeighbors(u, nbr_cnt);
        for (uint32_t j = 0; j < nbr_cnt; ++j) {
            uint32_t v = nbrs[j];
            uint32_t src = u;
            uint32_t dst = v;
            uint32_t kp = 0;
            if (src > dst) {
                std::swap(src, dst);
                kp = 1;
            }

            projection_operator->execute(&storage->edge_relations_[src][dst], kp, storage->candidate_sets_[u], storage->num_candidates_[u]);
            break;
        }
    }

    delete projection_operator;

    uint32_t n = query_graph->getVerticesCount();
    for (uint32_t i = 1; i < n; ++i) {
        uint32_t u = order[i];
        for (uint32_t j = 0; j < i; ++j) {
            uint32_t bn = order[j];
            if (query_graph->checkEdgeExistence(bn, u)) {
                if (i < core_vertices_cnt) {
                    convert_to_encoded_relation(storage, bn, u);
                }
                else {
                    convert_to_hash_relation(storage, bn, u);
                }
            }
        }
    }
}

void EvaluateQuery::convert_to_encoded_relation(catalog *storage, uint32_t u, uint32_t v) {
    uint32_t src = std::min(u, v);
    uint32_t dst = std::max(u, v);
    edge_relation& target_edge_relation = storage->edge_relations_[src][dst];
    edge* edges = target_edge_relation.edges_;
    uint32_t edge_size = target_edge_relation.size_;
    auto max_vertex_id = storage->data_graph_->getVerticesCount();
    assert(edge_size > 0);

    auto buffer = new uint32_t[max_vertex_id];
    memset(buffer, 0, sizeof(uint32_t)*max_vertex_id);

    uint32_t v_candidates_cnt = storage->get_num_candidates(v);
    uint32_t* v_candidates = storage->get_candidates(v);

    for (uint32_t i = 0; i < v_candidates_cnt; ++i) {
        uint32_t v_candidate = v_candidates[i];
        buffer[v_candidate] = i + 1;
    }

    uint32_t u_p = 0;
    uint32_t v_p = 1;
    if (u > v) {
        // Sort R(v, u) by u.
        std::sort(edges, edges + edge_size, [](edge& l, edge& r) -> bool {
            if (l.vertices_[1] == r.vertices_[1])
                return l.vertices_[0] < r.vertices_[0];
            return l.vertices_[1] < r.vertices_[1];
        });
        u_p = 1;
        v_p = 0;
    }

    encoded_trie_relation& target_encoded_trie_relation = storage->encoded_trie_relations_[u][v];
    uint32_t edge_cnt = edge_size;
    uint32_t u_candidates_cnt = storage->get_num_candidates(u);
    uint32_t* u_candidates = storage->get_candidates(u);
    target_encoded_trie_relation.size_ = u_candidates_cnt;
    target_encoded_trie_relation.offset_ = new uint32_t[u_candidates_cnt + 1];
    target_encoded_trie_relation.children_ = new uint32_t[edge_size];

    uint32_t offset = 0;
    uint32_t edge_index = 0;

    for (uint32_t i = 0; i < u_candidates_cnt; ++i) {
        uint32_t u_candidate = u_candidates[i];
        target_encoded_trie_relation.offset_[i] = offset;
        uint32_t local_degree = 0;
        while (edge_index < edge_cnt) {
            uint32_t u0 = edges[edge_index].vertices_[u_p];
            uint32_t v0 = edges[edge_index].vertices_[v_p];
            if (u0 == u_candidate) {
                if (buffer[v0] > 0) {
                    target_encoded_trie_relation.children_[offset + local_degree] = buffer[v0] - 1;
                    local_degree += 1;
                }
            }
            else if (u0 > u_candidate) {
                break;
            }

            edge_index += 1;
        }

        offset += local_degree;

        if (local_degree > target_encoded_trie_relation.max_degree_) {
            target_encoded_trie_relation.max_degree_ = local_degree;
        }
    }

    target_encoded_trie_relation.offset_[u_candidates_cnt] = offset;

    for (uint32_t i = 0; i < v_candidates_cnt; ++i) {
        uint32_t v_candidate = v_candidates[i];
        buffer[v_candidate] = 0;
    }
}

void EvaluateQuery::convert_to_hash_relation(catalog *storage, uint32_t u, uint32_t v) {
    // We assume that the relation is ordered.
    uint32_t src = std::min(u, v);
    uint32_t dst = std::max(u, v);
    edge_relation& target_edge_relation = storage->edge_relations_[src][dst];
    hash_relation& target_hash_relation1 = storage->hash_relations_[u][v];
    auto max_vertex_id = storage->data_graph_->getVerticesCount();

    edge* edges = target_edge_relation.edges_;
    uint32_t edge_size = target_edge_relation.size_;

    assert(edge_size > 0);

    uint32_t u_key = 0;
    uint32_t v_key = 1;

    if (src != u) {
        std::swap(u_key, v_key);
        // Sort the target edge relation.
        std::sort(edges, edges + edge_size, [](edge& l, edge& r)-> bool {
            if (l.vertices_[1] == r.vertices_[1])
                return l.vertices_[0] < r.vertices_[0];
            return l.vertices_[1] < r.vertices_[1];
        });
    }

    target_hash_relation1.children_ = new uint32_t[edge_size];

    uint32_t offset = 0;
    uint32_t local_degree = 0;
    uint32_t prev_u0 = max_vertex_id + 1;

    for (uint32_t i = 0; i < edge_size; ++i) {
        uint32_t u0 = edges[i].vertices_[u_key];
        uint32_t u1 = edges[i].vertices_[v_key];
        if (u0 != prev_u0 ) {
            if (prev_u0 != max_vertex_id + 1)
                target_hash_relation1.trie_->emplace(prev_u0, std::make_pair(local_degree, offset));

            offset += local_degree;

            if (local_degree > target_hash_relation1.max_degree_) {
                target_hash_relation1.max_degree_ = local_degree;
            }

            local_degree = 0;
            prev_u0 = u0;
        }

        target_hash_relation1.children_[offset + local_degree] = u1;
        local_degree += 1;
    }

    target_hash_relation1.cardinality_ = edge_size;
    target_hash_relation1.trie_->emplace(prev_u0, std::make_pair(local_degree, offset));
    if (local_degree > target_hash_relation1.max_degree_) {
        target_hash_relation1.max_degree_ = local_degree;
    }
}

void EvaluateQuery::convert_encoded_relation_to_sparse_bitmap(catalog *storage, ui*order) {
    uint32_t core_vertices_cnt = storage->query_graph_->get2CoreSize();

    for (uint32_t i = 1; i < core_vertices_cnt; ++i) {
        uint32_t u = order[i];

        for (uint32_t j = 0; j < i; ++j) {
            uint32_t bn = order[j];
            if (storage->query_graph_->checkEdgeExistence(u, bn)) {
                storage->bsr_relations_[bn][u].load(storage->encoded_trie_relations_[bn][u].get_size(),
                                                   storage->encoded_trie_relations_[bn][u].offset_,
                                                   storage->encoded_trie_relations_[bn][u].offset_,
                                                   storage->encoded_trie_relations_[bn][u].children_,
                                                   storage->max_num_candidates_per_vertex_, true);
            }
        }
    }
}

void EvaluateQuery::ComputeValidCans(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count,
                      ui**valid_cans, ui*valid_cans_count, ui* embedding, VertexID u, bool* visited_u, bool * visited_v) {
    ui unbrs_count;
    auto unbrs = query_graph->getVertexNeighbors(u, unbrs_count);
#ifdef ELABELED_GRAPH
    auto elabels = query_graph->getVertexEdgeLabels(u, unbrs_count);
#endif
    valid_cans_count[u] = 0;
    for (ui i = 0; i < candidates_count[u]; i++) {
        VertexID v = candidates[u][i];
        if (visited_v[v]) continue;
        bool flag = true;
        for (ui j = 0; j < unbrs_count; j++) {
            if (visited_u[unbrs[j]] == true
#ifdef ELABELED_GRAPH
                && !data_graph->checkEdgeExistence(v, embedding[unbrs[j]], elabels[j])) {
#else
                && !data_graph->checkEdgeExistence(v, embedding[unbrs[j]])) {
#endif
                flag = false;
                break;
            }
        }
        if (flag == true) {
            valid_cans[u][valid_cans_count[u]++] = v;
        }
    }
}

bool
EvaluateQuery::KSSEngine(const Graph *query_graph, const Graph *data_graph, Edges***edge_matrix,
                               ui **candidates, ui *candidates_count, ui *order, uint64_t output_limit_num,
                               uint64_t &call_cnt, mpz_t embedding_cnt, int64_t& time_limit) {
    ui cur_depth = 0;
    ui q_num = query_graph->getVerticesCount();
    ui d_num = data_graph->getVerticesCount();
    // 根据节点的度信息将order中的点，分为 kernel 和 shell 部分
    ui kernel_num = 0, shell_num = 0;
    ui* kernel = new ui[q_num];
    ui* shell = new ui[q_num];
    // kernel(true) or shell(false)
    bool* kos = new bool[q_num];
    memset(kos, 0, sizeof(bool)*q_num);
    // 划分shell & kernel: 沿着order做，每次判断下一个点是不是shell就好了
    ui * degree = new ui[q_num];
    for (ui i = 0; i < q_num; i++) {
        degree[i] = query_graph->getVertexDegree(i);
    }
    for (ui i = 0; i < q_num; i++) {
        VertexID u = order[i];
        if (degree[u] == 0) {
            shell[shell_num++] = u;
            continue;
        }
        kernel[kernel_num++] = u;
        kos[u] = true;
        ui nbr_num = 0;
        const ui* nbrs = query_graph->getVertexNeighbors(u, nbr_num);
        for (ui j = 0; j < nbr_num; j++) {
            degree[nbrs[j]]--;  // degree[nbrs[i]]--; 原来这里写错了...
        }
    }

    // 记录每一个shell相连的kernel点数量，每往下找一层，就检查一遍shell
    // 如果某个点的shell为0，计算这个shell的cans，如果cans为空，就回溯
    // 按照这样的算法，kernel计算完之后，shell也就算完了，直接算结果
    ui* shell2kernel = new ui[q_num];
    memset(shell2kernel, 0, sizeof(ui)*q_num);
    for (ui i = 0; i < shell_num; i++) {
        VertexID v = shell[i];
        ui nbr_num = 0;
        const ui* nbrs = query_graph->getVertexNeighbors(v, nbr_num);
        for (ui j = 0; j < nbr_num; j++) {
            if (kos[nbrs[j]] == true) {
                shell2kernel[v]++;
            }
        }
    }

    // allocate memory for auxiliary vars
    ui *idx = new ui [kernel_num];                   // depth as idx
    ui *idx_count = new ui [kernel_num];               // depth as idx
    ui *embedding = new ui [q_num];     // vid as idx
    ui **valid_cans = new ui* [q_num];  // vid as idx
    for (ui i = 0; i < q_num; i++) {
        valid_cans[i] = new ui [candidates_count[i]];// idx as idx
    }
    ui *valid_cans_cnt = new ui [q_num];// vid as idx
    bool *visited_v = new bool [d_num];  // vid as idx
    memset(visited_v, 0, sizeof(bool)*d_num);
    bool *visited_u = new bool[q_num];  // vid as idx
    memset(visited_u, 0, sizeof(bool)*q_num);

    bool overtime = false;
    VertexID start_vertex = kernel[0];
    visited_u[start_vertex] = true;
    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_cans[start_vertex][i] = candidates[start_vertex][i];
    }

    std::vector<VertexID> update; // (zhijie) 实际上, 对于固定的order, 不需要实时计算update

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = kernel[cur_depth];
            VertexID v = valid_cans[u][idx[cur_depth]];

            embedding[u] = v;
            visited_v[v] = true;
            idx[cur_depth] += 1;

            // 已经被计算的shell的cans检查
            // TODO: 要不要加?

            if (TimeOp::getClockNan() >= time_limit) {
                overtime = true;
                goto EXIT;
            }

            // 更新并检查shell2kernel的信息
            updateShell2Kernel(query_graph, u, shell2kernel, kos, update);

            if (cur_depth == kernel_num - 1) {
                for(ui i = 0; i < shell_num; i++) {
                    VertexID u_shell = shell[i];
                    ComputeValidCans(data_graph, query_graph, candidates, candidates_count, valid_cans,
                                 valid_cans_cnt, embedding, u_shell, visited_u, visited_v);
                }
                // 根据计算出的shell的cans信息，更新embedding_cnt
                if (computeKSSEmbeddingNaive(shell_num, shell, valid_cans, valid_cans_cnt, visited_v, embedding_cnt, time_limit) == true) {
                    overtime = true;
                    goto EXIT;
                }
                // embedding_cnt += computeKSSEmbeddingOpt1(shell_num, shell, valid_cans, valid_cans_cnt);
                // embedding_cnt += computeKSSEmbeddingOpt2(shell_num, shell, valid_cans, valid_cans_cnt);
                visited_v[v] = false;
                if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
                    goto EXIT;
                }
                // 最后一个节点匹配结束，恢复shell2kernel信息
                restoreShell2Kernel(query_graph, u, shell2kernel, kos);
            } else {
                cur_depth += 1;
                VertexID next_u = kernel[cur_depth];
                idx[cur_depth] = 0;
                call_cnt += 1;
                ComputeValidCans(data_graph, query_graph, candidates, candidates_count, valid_cans,
                                 valid_cans_cnt, embedding, next_u, visited_u, visited_v);
                
                visited_u[next_u] = true;
                idx_count[cur_depth] = valid_cans_cnt[next_u];
            }
        }

        // 回溯部分
        cur_depth -= 1;
        if (cur_depth == (ui)-1)
            break;
        VertexID last_u = kernel[cur_depth + 1];
        visited_v[embedding[kernel[cur_depth]]] = false;
        visited_u[last_u] = false;
        // 恢复shell2kernel信息
        restoreShell2Kernel(query_graph, kernel[cur_depth], shell2kernel, kos);
    }

    // Release the memory
    EXIT:
    delete[] kernel;
    delete[] shell;
    delete[] kos;
    delete[] degree;
    delete[] shell2kernel;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    for (ui i = 0; i < q_num; i++) {
        delete[] valid_cans[i];
    }
    delete[] valid_cans;
    delete[] valid_cans_cnt;
    delete[] visited_u;
    delete[] visited_v;
    return overtime;
}

void // 更新shell2kernel信息 // (zhijie) 能不能传引用呀 (以避免重复)
EvaluateQuery::updateShell2Kernel(const Graph *query_graph, VertexID u, ui* shell2kernel, bool* kos, std::vector<VertexID> & update) {
    ui nbr_num = 0;
    const ui* nbrs = query_graph->getVertexNeighbors(u, nbr_num);
    update.clear();
    for (ui i = 0; i < nbr_num; i++) {
        VertexID nbr = nbrs[i];
        if (kos[nbr] == false) {
            shell2kernel[nbr]--;
            if(shell2kernel[nbr] == 0)
                update.emplace_back(nbr);
        }
    }
}

void  // 恢复shell2kernel信息
EvaluateQuery::restoreShell2Kernel(const Graph *query_graph, VertexID u, ui* shell2kernel, bool* kos) {
    ui nbr_num = 0;
    const ui*nbrs = query_graph->getVertexNeighbors(u, nbr_num);
    for (ui i = 0; i < nbr_num; i++) {
        VertexID nbr = nbrs[i];
        if (kos[nbr] == false) {
            shell2kernel[nbr]++;
        }
    }
}

// Naive 方法
bool 
EvaluateQuery::computeKSSEmbeddingNaive(ui shell_num, ui* shell, ui** valid_cans, ui* valid_cans_count, bool * visited_v,
                                        mpz_t embedding_cnt, int64_t& time_limit) {
    
    return computeKSSEmbeddingNaiveImpl(0, shell_num, shell, valid_cans, valid_cans_count, visited_v, embedding_cnt, time_limit);
}

bool 
EvaluateQuery::computeKSSEmbeddingNaiveImpl(ui depth, ui shell_num, ui* shell, ui** valid_cans, ui* valid_cans_count, bool * visited_v,
                                            mpz_t embedding_cnt, int64_t& time_limit) {

    if (TimeOp::getClockNan() >= time_limit) {
        return true;
    }
    VertexID u_shell = shell[depth];

    // 对于叶子节点, 避免展开
    if (depth == shell_num - 1) {
        // 使用 visited_u, 无法用 valid_cans_count[u_shell] - dup 的方式对叶节点计算
        for(ui i = 0; i < valid_cans_count[u_shell]; i++) {
            VertexID v_id = valid_cans[u_shell][i];
            if (!visited_v[v_id]) mpz_add_ui(embedding_cnt, embedding_cnt, 1);
        }
    }
    // 对于剩余的节点，不加检查的进行展开 [与树状节点的区别是, 甚至不需要对邻居进行获取]
    else {
        for(ui i = 0; i < valid_cans_count[u_shell]; i++) {
            VertexID v_id = valid_cans[u_shell][i];
            if(!visited_v[v_id]){
                visited_v[v_id] = true;
                if (computeKSSEmbeddingNaiveImpl(depth+1, shell_num, shell, valid_cans, valid_cans_count, visited_v,
                                                 embedding_cnt, time_limit) == true)
                    return true;
                visited_v[v_id] = false;
            }
        }
    }
    return false;
}

#if BSX_SIM_THRESHOLD != 100
SimInfo* BatchInfo::sim_info = nullptr;
#endif

/**
 * use bsx method
 * not support edge label(so far)
*/
bool
EvaluateQuery::BSXEngine(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix,
                          ui **candidates, ui *candidates_count, uint64_t output_limit_num,
                          uint64_t &call_count, mpz_t embedding_cnt, int64_t& time_limit) {
#ifdef ELABELED_GRAPH
    return 0;
#endif
    ui q_num = query_graph->getVerticesCount();
    ui max_candidates_num = candidates_count[0];
    for (ui i = 1; i < q_num; ++i) {
        if (candidates_count[i] > max_candidates_num) {
            max_candidates_num = candidates_count[i];
        }
    }
    ui* order = nullptr;
    // separate leaf and trunk vertices(min_vertex_cover)
    ui num_cover = 0;
    bsxMaxCoverOrder(query_graph, order, num_cover, candidates_count);
    auto num_indep = q_num - num_cover;
    const VertexID* indep_nodes =  order + num_cover;

    // structure used to store history intersection info
    // std::deque<IntersectCache> cachedIntersect;   // max size is the height of the tree
    // construct index structure, contains the history info
    // new index structure, update in time, 24-3-7
    BSXIndex index(query_graph, data_graph, edge_matrix, candidates, candidates_count, num_cover);
    auto& batch_info = index.batch_info;

    // auxiliary data structure
    bool overtime = false;
    auto& visited_u = index.visited_u;
    auto& u2v = index.embedding->u2v;
    auto& depth2u = index.embedding->depth2u;
    mpz_init_set_ui(embedding_cnt, 0);
    ui cur_depth = 0;
    VertexID start_vertex = order[cur_depth];
    depth2u.emplace_back(start_vertex);
    visited_u[start_vertex] = true;
    auto& level_embeddings = index.level_embeddings_;

#if BSX_SIM_THRESHOLD != 100
    BatchInfo::sim_info = new SimInfo();
    auto& enable_sb = BatchInfo::sim_info->enable_;
    auto& sim_batches = BatchInfo::sim_info->sim_batches_;
    auto& sb_offset = BatchInfo::sim_info->offset_;
    auto& sb_num = BatchInfo::sim_info->num_;
    auto& sb_idxs = BatchInfo::sim_info->idxs_;
    auto& sb_valid = BatchInfo::sim_info->valid_;
    auto& sb_max_diff = BatchInfo::sim_info->max_diff_;
    auto& sb_diff_unbrs = BatchInfo::sim_info->diff_unbrs_;
    auto& sb_diff_unbrs_cnt = BatchInfo::sim_info->diff_unbrs_cnt_;
    auto& sb_same_or_diff = BatchInfo::sim_info->same_or_diff_;
    auto& sb_judge_nbrs = BatchInfo::sim_info->judge_nbrs_;
    auto& batches_processed = BatchInfo::sim_info->batches_processed_;
    if ((query_graph->getGraphMaxDegree() - 1) * (1 - BSX_SIM_THRESHOLD * 0.01) >= 1) {
        // allowed max different u_nbrs is less than 1
        enable_sb = true;
    }
    if (enable_sb) {
        BatchInfo::sim_info->init(query_graph, max_candidates_num);
    }
#endif

    // init info of start vertex
    batch_info[start_vertex].add();
    bsxComEqBatch(index, start_vertex);
    index.valid_cans_[start_vertex].push(new VertexID[batch_info[start_vertex].maxCnt_.top()]);
    index.valid_cnt_[start_vertex].push(0);

#if BSX_SIM_THRESHOLD != 100
    if (enable_sb) {
        if (sb_max_diff[start_vertex]) bsxComSimBatch(index, start_vertex);
    }
#endif

#ifdef ANALYZE_DUPLICATE
    auto g_name = query_graph->g_name;
    size_t last_slash_pos = g_name.find_last_of('/');
    if (last_slash_pos != std::string::npos)
        g_name = g_name.substr(last_slash_pos + 1);
    size_t last_dot_pos = g_name.find_last_of('.');
    if (last_dot_pos != std::string::npos)
        g_name = g_name.substr(0, last_dot_pos);
    g_name = "./" + g_name;
    int status = mkdir(g_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if (status != 0) {
        std::cout << g_name << ": Failed to create folder.\n";
        exit(-1);
    }

    std::vector<std::ofstream> out_files;
    out_files.resize(q_num);
    for (ui i = 1; i < q_num; i++) {
        out_files[i-1].open(g_name + "/" + std::to_string(i-1) + ".txt");
    }
#endif

    while (true) {
        while (batch_info[depth2u[cur_depth]].idx_.top() < batch_info[depth2u[cur_depth]].num_.top()) {
            VertexID u = depth2u[cur_depth];
#if BSX_SIM_THRESHOLD != 100
            // do not judge sb_valid(just store current state)
            // do not process last cover_vertex
            if (enable_sb && sb_max_diff[u] && cur_depth < num_cover - 1) {
                while(batches_processed[u][batch_info[depth2u[cur_depth]].idx_.top()]
                      && batch_info[depth2u[cur_depth]].idx_.top() < batch_info[depth2u[cur_depth]].num_.top()) {
                    batch_info[depth2u[cur_depth]].idx_.top()++;
                }
                if (batch_info[depth2u[cur_depth]].idx_.top() >= batch_info[depth2u[cur_depth]].num_.top()) {
                    break;
                }
                batches_processed[u][batch_info[depth2u[cur_depth]].idx_.top()] = true;
            }
#endif
            ui cur_batch_cnt;
            VertexID* cur_batch = batch_info[u].cur_batch(cur_batch_cnt);
            if (TimeOp::getClockNan() >= time_limit) {
                overtime = true;
                goto EXIT;
            }

            // nxt batch
            batch_info[depth2u[cur_depth]].idx_.top()++;

            auto& cur_cans_cnt = index.valid_cnt_[u].top();
            auto& cur_cans = index.valid_cans_[u].top();
            std::copy(cur_batch, cur_batch+cur_batch_cnt, cur_cans);
            cur_cans_cnt = cur_batch_cnt;

            VertexID failed_u = bsxRefine(index, u);
            if (failed_u != (ui)-1) {  // no valid cans for next depth
#if BSX_SIM_THRESHOLD != 100
                if (enable_sb) {
                    // process the diff_nbrs of fail_u
                    for (ui i = 0; i < sb_diff_unbrs_cnt[failed_u]; i++) {
                        auto& sb_diff_unbr = sb_diff_unbrs[failed_u][i];
                        sb_valid[sb_diff_unbr] = false;
                    }
                }
#endif
                continue;
            }

#ifdef ANALYZE_DUPLICATE
            for (ui i = 0; i < q_num; i++) {
                auto cur_u = order[i];
                if (!visited_u[cur_u]) {
                    for (ui j = 0; j < index.valid_cnt_[cur_u].top();j++) {
                        out_files[cur_depth] << index.valid_cans_[cur_u].top()[j] << " ";
                    }
                    out_files[cur_depth] << std::endl;
                }
            }
            out_files[cur_depth] << "------" << std::endl;
#endif

            u2v[u] = cur_batch[0];

            if (cur_depth >= num_cover - 1) {
#ifdef ANALYZE_DUPLICATE
                for (ui indep_idx = num_cover; indep_idx < q_num - 1; indep_idx++) {
                    for (ui i = num_cover; i < q_num; i++) {
                        auto cur_u = order[i];
                        for (ui j = 0; j < index.valid_cnt_[cur_u].top();j++) {
                            out_files[indep_idx] << index.valid_cans_[cur_u].top()[j] << " ";
                        }
                        out_files[indep_idx] << std::endl;
                    }
                    out_files[indep_idx] << "------" << std::endl;
                }
#endif
                // enumerate results on indep nodes, process ancestors' ves by the way
                bsxGenResult(num_indep, indep_nodes, index);
                mpz_add(embedding_cnt, embedding_cnt, level_embeddings);
                if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
                    goto EXIT;
                }
                // next batch
                bsxDeRefine(index);
#if BSX_SIM_THRESHOLD != 100
                // do not process vertex conlicts
                if (enable_sb) {
                    memset(sb_valid, false, sizeof(bool)*q_num);
                }
#endif
            } else {
                cur_depth++;
                VertexID cur_u = bsxGenNxtU(index, order, cur_depth, num_cover);
                if (cur_u == (VertexID)-1) cur_u = order[cur_depth];
                depth2u.emplace_back(cur_u);
                call_count++;
                // construct nbrs&seperate batches, and then refinement
                batch_info[cur_u].add();
                bsxComEqBatch(index, cur_u);
                index.valid_cans_[cur_u].push(new VertexID[batch_info[cur_u].maxCnt_.top()]);
                index.valid_cnt_[cur_u].push(0);
                visited_u[cur_u] = true;
#if BSX_SIM_THRESHOLD != 100
                if (enable_sb) {
                    if (sb_max_diff[cur_u] && cur_depth < num_cover - 1)
                        bsxComSimBatch(index, cur_u);
                    // process the diff_nbrs of cur_u
                    for (ui i = 0; i < sb_diff_unbrs_cnt[cur_u]; i++) {
                        auto& sb_diff_unbr = sb_diff_unbrs[cur_u][i];
                        sb_valid[sb_diff_unbr] = false;
                    }
                }
#endif
            }
        }

        // backtracking
        cur_depth -= 1;
        if (cur_depth == ui(-1))
            break;
        VertexID last_u = depth2u[cur_depth+1];
        VertexID cur_u = depth2u[cur_depth];
        depth2u.resize(cur_depth+1);
        visited_u[last_u] = false;

#if BSX_SIM_THRESHOLD != 100  // clear the diff nbrs
        if (enable_sb && sb_max_diff[last_u] && cur_depth + 1 < num_cover - 1) {
            auto last_b_num = batch_info[last_u].num_.top();
            if (last_b_num != 0) {
                ui unbrs_count;
                const ui *unbrs = index.q_graph_->getVertexNeighbors(last_u, unbrs_count);
                ui need_cnt = unbrs_count - sb_max_diff[last_u];
                for (ui i = 0; i < unbrs_count; i++) {
                    auto& unbr = unbrs[i];
                    if (visited_u[unbr]) continue;
                    if (!need_cnt) {
                        sb_diff_unbrs_cnt[unbr]--;
                    } else {
                        need_cnt--;
                    }
                }
            }
            if (sb_valid[cur_u]) {
                auto b_idx = batch_info[cur_u].idx_.top() - 1;
                auto sb_idx = sb_idxs[cur_u][b_idx] - 1;
                auto sb_cnt = sb_offset[cur_u][sb_idx + 1] - sb_offset[cur_u][sb_idx];
                for (ui i = 0; i < sb_cnt; i++) {
                    auto& batch_idx = sim_batches[cur_u][sb_offset[cur_u][sb_idx]+i];
                    batches_processed[cur_u][batch_idx] = true;
                }
            }
        }
#endif
        batch_info[last_u].pop();
        delete[] index.valid_cans_[last_u].top();
        index.valid_cans_[last_u].pop();
        index.valid_cnt_[last_u].pop();
        bsxDeRefine(index);
    }

    // Release the buffer.
    EXIT:

#if BSX_SIM_THRESHOLD != 100
    delete BatchInfo::sim_info;
#endif

#ifdef ANALYZE_DUPLICATE
    for (ui i = 1; i < q_num; i++) {
        out_files[i-1].close();
    }
#endif

    return overtime;
}

/**
 * generate order for BSXEngine
 * 1. static: seperate vertices into cover&independent by MaxCover
 * 2. dynamic: sort each kind by #cans(asc), #degree(des), #id(asc)
 * compute static order by max indep cover, dynamic order computed along matching
*/
void
EvaluateQuery::bsxMaxCoverOrder(const Graph *graph, ui *&order, ui& num_cover, ui *candidates_count) {
    auto q_num = graph->getVerticesCount();
    if (order == nullptr) {
        order = new ui[q_num];
    }
    IndepSet indepSet(graph);
    auto indep = indepSet.linearTime();
    num_cover = q_num - indep.second;
    // complete order
    ui num_indep = 0;
    ui num_other = 0;
    for (ui i = 0; i < q_num; i++) {
        if (indep.first[i]) {
            order[num_cover + num_indep++] = i;
        } else {
            order[num_other++] = i;
        }
    }
    // compute the first u
    ui first_u = order[0];
    ui first_idx = 0;
    for (ui i = 1; i < num_cover; i++) {
        ui cur_u = order[i];
        if (candidates_count[cur_u] < candidates_count[first_u]
            || (candidates_count[cur_u] == candidates_count[first_u]
                && graph->getVertexDegree(cur_u) > graph->getVertexDegree(first_u))) {
            first_idx = i;
            first_u = cur_u;
        }
    }
    order[first_idx] = order[0];
    order[0] = first_u;

    // dynamic->static order
    // std::sort(order, order+num_cover,[candidates_count, graph](VertexID l, VertexID r){
    //     if (candidates_count[l] == candidates_count[r]) {
    //         if (graph->getVertexDegree(l) == graph->getVertexDegree(r)) {
    //             return l < r;  // id(asc)
    //         }
    //         return graph->getVertexDegree(l) > graph->getVertexDegree(r);  // degree(desc)
    //     }
    //     return candidates_count[l] < candidates_count[r];  // cans(asc)
    // });
    delete[] indep.first;
}

/**Reverse op of BSXRefine
 * pop the valid_cans of influenced_u
 * process oneCansV from refinement
*/
void
EvaluateQuery::bsxDeRefine(BSXIndex& index) {
    auto q_num = index.q_graph_->getVerticesCount();
    bool* nbr_updated = new bool[q_num];
    memset(nbr_updated, false, sizeof(bool)*q_num);

    // process first u in influenced_u seperately, valid_cans of first_influenced_u comes from
    //   its batch_nodes, and couldn't be deleted
    auto first_u = index.influenced_u_.top()[0];
    nbr_updated[first_u] = true;
    // remove index_, added at refinement
    ui nbrs_cnt;
    auto nbrs = index.q_graph_->getVertexNeighbors(first_u, nbrs_cnt);
    for (ui i = 0; i < nbrs_cnt; i++) {
        auto& nbr = nbrs[i];
        delete index.index_[nbr][first_u].top();
        delete index.index_[first_u][nbr].top();
        index.index_[nbr][first_u].pop();
        index.index_[first_u][nbr].pop();
    }
    // recover index_cans_, changed at refinement
    auto tmp_cans = index.valid_cans_[first_u].top();
    index.valid_cans_[first_u].pop();
    index.valid_cnt_[first_u].pop();
    index.index_cans_[first_u] = index.valid_cans_[first_u].top();
    index.index_cnt_[first_u] = index.valid_cnt_[first_u].top();
    index.valid_cans_[first_u].push(tmp_cans);
    index.valid_cnt_[first_u].push(0);

    for (ui i = 1; i < index.influenced_u_.top().size(); i++) {
        auto influenced_u = index.influenced_u_.top()[i];
        // not delete first influenced_u, because its valid_cans is not constructed in refine
        delete[] index.valid_cans_[influenced_u].top();
        index.valid_cans_[influenced_u].pop();
        index.valid_cnt_[influenced_u].pop();
        index.index_cans_[influenced_u] = index.valid_cans_[influenced_u].top();
        index.index_cnt_[influenced_u] = index.valid_cnt_[influenced_u].top();

        // remove index_
        ui nbrs_cnt;
        auto nbrs = index.q_graph_->getVertexNeighbors(influenced_u, nbrs_cnt);
        for (ui i = 0; i < nbrs_cnt; i++) {
            auto& nbr = nbrs[i];
            if (nbr_updated[nbr]) continue;
            delete index.index_[nbr][influenced_u].top();
            delete index.index_[influenced_u][nbr].top();
            index.index_[nbr][influenced_u].pop();
            index.index_[influenced_u][nbr].pop();
        }
        nbr_updated[influenced_u] = true;
    }
    index.influenced_u_.pop();
    delete[] nbr_updated;
}

/**generate next u
 * sort by:
 * 1.depth < num_cover: not_matched, #cans(asc), #degree(des), #id(arbitrary)
 * 2.depth >= num_cover: #cans!=1, #cans(asc), #degree(des), #id(arbitrary)
*/
VertexID
EvaluateQuery::bsxGenNxtU(BSXIndex& index, VertexID* order, ui depth, ui num_cover) {
    auto& valid_cnt = index.valid_cnt_;
    auto& graph = index.q_graph_;
    if (depth < num_cover) {
        std::sort(order+depth, order+num_cover, [valid_cnt, graph](VertexID a, VertexID b) {
            if (valid_cnt[a].top() == valid_cnt[b].top()) {
                return graph->getVertexDegree(a) > graph->getVertexDegree(b);
            }
            return valid_cnt[a].top() < valid_cnt[b].top();
        });
        return (VertexID)-1;
    } else {
        VertexID* tmp_nodes = new VertexID[num_cover];
        std::copy(order,order+num_cover, tmp_nodes);
        std::sort(tmp_nodes, tmp_nodes+num_cover, [valid_cnt, graph](ui a, ui b) {
            if (valid_cnt[a].top() != 1 && valid_cnt[b].top() != 1) {
                if (valid_cnt[a].top() == valid_cnt[b].top()) {
                    return graph->getVertexDegree(a) > graph->getVertexDegree(b);
                }
                return valid_cnt[a].top() < valid_cnt[b].top();
            }
            return valid_cnt[a].top() > valid_cnt[b].top();
        });
        VertexID node = tmp_nodes[0];
        delete[]tmp_nodes;
        return node;
    }
    return (VertexID)-1;
}

// check termination (each no-indep only one cans)
bool
EvaluateQuery::bsxCheckTermination(ui num, VertexID* indep, std::stack<ui>*valid_cnt) {
    while (num) if (valid_cnt[indep[--num]].top() != 1) return false;
    return true;
}

// compute valid_cans for all indep, detect conflict
bool
EvaluateQuery::bsxGenIndepValidCans(ui indep_num, const VertexID* indep, BSXIndex& index, std::vector<std::vector<VertexID>>& cans) {
    const VertexID** uu_nbrs = new const VertexID*[index.q_graph_->getVerticesCount()];
    ui* uu_nbrs_cnt = new  ui[index.q_graph_->getVerticesCount()];
    for (ui i = 0; i < indep_num; i++) {
        auto u = indep[i];
        ui nbrs_cnt = 0;
        auto nbrs = index.q_graph_->getVertexNeighbors(u, nbrs_cnt);
        for (ui j = 0; j < nbrs_cnt; j++) {
            auto& nbr = nbrs[j];
            auto& nbr_v = index.embedding->u2v[nbr];
            uu_nbrs[j] = index.getNeighbors(nbr, u, nbr_v, uu_nbrs_cnt[j]);
        }
        auto intersected = std::move(SetOp::intersectMultiple(uu_nbrs, uu_nbrs_cnt, nbrs_cnt));
        if (intersected.size() == 0) {
            delete[] uu_nbrs;
            delete[] uu_nbrs_cnt;
            return false;
        }
        cans[u] = std::move(intersected);
    }
    delete[] uu_nbrs;
    delete[] uu_nbrs_cnt;
    return true;
}

// generate equivalent batches
void
EvaluateQuery::bsxComEqBatch(BSXIndex& index, VertexID u) {
    auto& num_node = index.valid_cnt_[u].top();
    auto& batch_nodes = index.batch_info[u].nodes_.top();
    auto& offset = index.batch_info[u].offset_.top();
    auto& cnt = index.batch_info[u].cnt_.top();
    batch_nodes = new VertexID[num_node];
    offset = new VertexID[num_node];
    cnt = new VertexID[num_node];
    if (index.valid_cnt_[u].top() < 16) {  // if #nodes < 16, compute batch by comparing each other
        std::vector<ui> idxs;
        idxs.reserve(num_node);
        for (ui i = 0; i < num_node; i++) idxs.emplace_back(i);
        bsxComEqBatchDirect(index, u, idxs);
    } else {  // if #nodes >= 16, compute batch by nodeSimilarity
        std::vector<int64_t> similarity = std::move(NodeSim::nodeSim(index, u));
        std::vector<std::pair<int64_t, ui>> sim_sorted;
        sim_sorted.reserve(num_node);
        for (ui i = 0; i < similarity.size(); i++) {
            sim_sorted.emplace_back(similarity[i], i);
        }
        // sort smilarity, id (ensure the correct order)
        std::sort(sim_sorted.begin(), sim_sorted.end());
        ui batch_start = 0, batch_end = 0;
        std::vector<VertexID> batch_idxs;
        while(batch_end < num_node) {
            batch_idxs.clear();
            batch_idxs.emplace_back(sim_sorted[batch_start].second);
            while((++batch_end) < num_node && sim_sorted[batch_end].first == sim_sorted[batch_start].first) {
                batch_idxs.emplace_back(sim_sorted[batch_end].second);;
            }
            assert(batch_end-batch_start == batch_idxs.size());
            bsxComEqBatchDirect(index, u, batch_idxs);
            batch_start = batch_end;
        }
    }
}

// compute equ-batch on idxs, idxs indicate which nodes participate batch computation
void
EvaluateQuery::bsxComEqBatchDirect(BSXIndex& index, VertexID u, std::vector<ui>& idxs) {
    auto& nodes = index.valid_cans_[u].top();
    auto num_idxs = idxs.size();
    auto& batches = index.batch_info[u].nodes_.top();
#if BSX_SIM_THRESHOLD != 100
    auto& enable_sb = BatchInfo::sim_info->enable_;
#endif
    auto& offset = index.batch_info[u].offset_.top();
    auto& cnt = index.batch_info[u].cnt_.top();
    auto& num = index.batch_info[u].num_.top();
    auto& maxCnt = index.batch_info[u].maxCnt_.top();
    ui unbrs_count;
    const ui *unbrs = index.q_graph_->getVertexNeighbors(u, unbrs_count);
    // -1 -- error-node(offset[*] == 0) and should be deleted, 0 -- un-processed
    ui* batch_idx = new ui[num_idxs];  // indicate each node belonging to ?th batch
    memset(batch_idx, 0, sizeof(ui)*num_idxs);

    // delete all nodes that have no edge connection with its neighbors
    for (ui i = 0; i < num_idxs; i++) {
        auto idx = idxs[i];
        for (ui unbr_idx = 0; unbr_idx < unbrs_count; unbr_idx++) {
            auto& unbr = unbrs[unbr_idx];
            auto& edges = index.index_[u][unbr].top();
            if (edges->offset_[idx+1] - edges->offset_[idx] == 0) {
#if BSX_SIM_THRESHOLD != 100
                if (enable_sb) {
                    auto& sb_valid = BatchInfo::sim_info->valid_;
                    auto& sb_diff_unbrs = BatchInfo::sim_info->diff_unbrs_;
                    auto& sb_diff_unbrs_cnt = BatchInfo::sim_info->diff_unbrs_cnt_;
                    // process the diff_nbrs of fail_u
                    for (ui i = 0; i < sb_diff_unbrs_cnt[unbr]; i++) {
                        auto& sb_diff_unbr = sb_diff_unbrs[unbr][i];
                        sb_valid[sb_diff_unbr] = false;
                    }
                }
#endif
                batch_idx[i] = (ui)-1;
                break;
            }
        }
    }

    ui batch_cnt = 1;  // number from 1
    for (ui i = 0; i < num_idxs; i++) {
        auto& idx = idxs[i];
        if (batch_idx[i] != 0) continue;
        batch_idx[i] = batch_cnt;
        offset[num] = num == 0 ? 0 : offset[num-1]+cnt[num-1];  // set offset of cur_batch
        cnt[num] = 1;
        batches[offset[num]] = nodes[idx];
#if BSX_SIM_THRESHOLD != 100
        if (enable_sb) {
            auto& enable_u_sb = BatchInfo::sim_info->max_diff_[u];
            auto& batches_idxs = BatchInfo::sim_info->batches_idxs_[u];
            if (enable_u_sb) batches_idxs[offset[num]] = idx;
        }
#endif
        for (ui j = i+1; j < num_idxs; j++) {
            if (batch_idx[j] != 0) continue;
            auto ev_idx = idxs[j];
            bool equ = true;
            // compare the nbrs
            for (ui unbr_idx = 0; unbr_idx < unbrs_count; unbr_idx++) {
                VertexID unbr = unbrs[unbr_idx];
                auto& edges = index.index_[u][unbr].top();
                if (edges->offset_[idx+1]-edges->offset_[idx]
                    != edges->offset_[ev_idx+1] - edges->offset_[ev_idx]) {
                    equ = false;
                    break;
                }
                for (ui u2 = 0; u2 < edges->offset_[idx+1] - edges->offset_[idx]; u2++) {
                    if (edges->edge_[u2+edges->offset_[idx]] != edges->edge_[u2+edges->offset_[ev_idx]]) {
                        equ = false;
                        break;
                    }
                }
                if (equ == false) break;
            }
            if (equ == true) {
                batch_idx[j] = batch_cnt;
                batches[offset[num]+cnt[num]] = nodes[ev_idx];
#if BSX_SIM_THRESHOLD != 100
                if (enable_sb) {
                    auto& enable_u_sb = BatchInfo::sim_info->max_diff_[u];
                    auto& batches_idxs = BatchInfo::sim_info->batches_idxs_[u];
                    if (enable_u_sb) batches_idxs[offset[num]+cnt[num]] = ev_idx;
                }
#endif
                cnt[num]++;
            }
        }
        // process max_cnt
        if (maxCnt < cnt[num]) maxCnt = cnt[num];
        // no need for sort, because the idxs have been sorted and this fun. will not break it
        // std::sort(batch+offset[num], batch+(offset[num]+cnt[num]));
        num++;
        batch_cnt++;
    }
    delete[] batch_idx;
}

// compute sim_batches based on the equal batches
void
EvaluateQuery::bsxComSimBatch(BSXIndex& index, VertexID u) {
#if BSX_SIM_THRESHOLD != 100
    static auto q_num = index.q_graph_->getVerticesCount();
    static ui judge_nbrs_cnt;
    auto& max_diff = BatchInfo::sim_info->max_diff_[u];
    auto& sb_valid = BatchInfo::sim_info->valid_;

    auto& b_idxs = BatchInfo::sim_info->batches_idxs_[u];
    auto b_num = index.batch_info[u].num_.top();
    auto b_fst = index.batch_info[u].offset_.top();  // use the 1st as the representation

    if (max_diff == 0 || b_num == 0) return;
    sb_valid[u] = true;

    auto& sim_batches = BatchInfo::sim_info->sim_batches_[u];
    auto& offset = BatchInfo::sim_info->offset_[u];
    auto& num = BatchInfo::sim_info->num_[u];
    auto& sb_idxs = BatchInfo::sim_info->idxs_[u];
    auto& diff_unbrs = BatchInfo::sim_info->diff_unbrs_;
    auto& diff_unbrs_cnt = BatchInfo::sim_info->diff_unbrs_cnt_;
    auto& visited_u = index.visited_u;
    auto& same_or_diff = BatchInfo::sim_info->same_or_diff_[u];
    memset(same_or_diff, false, sizeof(bool)*q_num);
    auto& judge_nbrs = BatchInfo::sim_info->judge_nbrs_;
    judge_nbrs_cnt = 0;
    auto& batches_processed = BatchInfo::sim_info->batches_processed_[u];
    memset(batches_processed, false, sizeof(bool)*b_num);
    // compute the nbrs that should be the same(with the same cans in these nbrs)
    ui matched_nbrs = 0;
    ui unbrs_count;
    const ui *unbrs = index.q_graph_->getVertexNeighbors(u, unbrs_count);
    ui need_cnt = unbrs_count - max_diff;
    for (ui i = 0; i < unbrs_count; i++) {
        auto& unbr = unbrs[i];
        if (visited_u[unbr]) {
            matched_nbrs++;
            continue;
        }
        if (!need_cnt) {
            same_or_diff[unbr] = true;
            diff_unbrs[unbr][diff_unbrs_cnt[unbr]++] = u;
        } else {
            need_cnt--;
            judge_nbrs[judge_nbrs_cnt++] = unbr;
        }
    }
    if (judge_nbrs_cnt + matched_nbrs == unbrs_count) {
        sb_valid[u] = false;
        return;
    }

    // 0 -- un-processed
    memset(sb_idxs, 0, sizeof(ui)*b_num);

    num = 1;
    offset[0] = 0;
    for (ui i = 0; i < b_num; i++) {
        auto& b_idx = b_idxs[b_fst[i]];
        if (sb_idxs[i] != 0) continue;
        sb_idxs[i] = num;
        offset[num] = offset[num-1];  // set offset of cur_batch
        sim_batches[offset[num]++] = i;
        for (ui j = i+1; j < b_num; j++) {
            if (sb_idxs[j] != 0) continue;
            auto sim_idx = b_idxs[b_fst[j]];
            bool sim = true;
            // compare the nbrs
            for (ui unbr_idx = 0; unbr_idx < judge_nbrs_cnt; unbr_idx++) {
                VertexID unbr = judge_nbrs[unbr_idx];
                auto& edges = index.index_[u][unbr].top();
                if (edges->offset_[b_idx+1]-edges->offset_[b_idx]
                    != edges->offset_[sim_idx+1] - edges->offset_[sim_idx]) {
                    sim = false;
                    break;
                }
                for (ui u2 = 0; u2 < edges->offset_[b_idx+1] - edges->offset_[b_idx]; u2++) {
                    if (edges->edge_[u2+edges->offset_[b_idx]] != edges->edge_[u2+edges->offset_[sim_idx]]) {
                        sim = false;
                        break;
                    }
                }
                if (sim == false) break;
            }
            if (sim == true) {
                sb_idxs[j] = num;
                sim_batches[offset[num]++] = j;
            }
        }
        num++;
    }
#endif
    return;
}

// equ-batch refine, just process first v of valid_cans, because of they are equ
// if success, return -1, else return failed uId
ui
EvaluateQuery::bsxRefine(BSXIndex& index, VertexID u) {
    ui q_num = index.q_graph_->getVerticesCount();
    auto v = index.valid_cans_[u].top()[0];
    bool* influenced = new bool[q_num];
    memset(influenced, false, sizeof(bool)*q_num);
    std::vector<VertexID> cur_inf;
    // if a node is influenced for the first time, just use the generated result
    // do not need intersection operation(op) with old valid_cans
    std::pair<const VertexID*, ui>* influenced_cans = new std::pair<const VertexID*, ui>[q_num];
    ui unbrs_cnt;
    auto unbrs = index.q_graph_->getVertexNeighbors(u, unbrs_cnt);
    ui returned_value = (ui)-1;
    for (ui i = 0; i < unbrs_cnt; i++) {
        auto& unbr = unbrs[i];
        if (index.visited_u[unbr]) continue;
        // old valid_cans of unbr
        auto& unbr_valid_cans = index.valid_cans_[unbr].top();
        auto& unbr_valid_cnt = index.valid_cnt_[unbr].top();
        ui vnbr_cnt;
        auto vnbrs = index.getNeighbors(u, unbr, v, vnbr_cnt);
        if (vnbr_cnt == 0) {
            returned_value = unbr;
            goto bsxRefine_EXIT;
        }
        // // vnbrs must be included in unbr_valid_cans
        // assert(SetOp::setInclude(vnbrs, vnbr_cnt, unbr_valid_cans, unbr_valid_cnt));

        if (unbr_valid_cnt == vnbr_cnt) continue;  // means no changes

        // stack a valid_cans&valid_cnt
        influenced[unbr] = true;
        influenced_cans[unbr] = std::make_pair(vnbrs, vnbr_cnt);
    }

    // write influenced valid_cans to index, then update index
    for (ui unbr = 0; unbr < q_num; unbr++) {
        if (influenced[unbr]) {
            auto& vnbr_cnt = influenced_cans[unbr].second;
            auto& vnbrs = influenced_cans[unbr].first;
            auto new_valid_cans = new VertexID[vnbr_cnt];
            std::copy(vnbrs, vnbrs+vnbr_cnt, new_valid_cans);
            index.valid_cans_[unbr].push(new_valid_cans);
            index.valid_cnt_[unbr].push(vnbr_cnt);
        }
    }
    influenced[u] = true;

    // implement index.updateOneSide(influenced) later
    index.updateIndex(influenced);
    // if (index.updateIndex(influenced) == (ui)-1) {
    //     influenced[u] = false;
    //     for (ui i = 0; i < q_num; i++) {
    //         if (influenced[i]) {
    //             delete[] index.valid_cans_[i].top();
    //             index.valid_cans_[i].pop();
    //             index.valid_cnt_[i].pop();
    //         }
    //     }
    //     returned_value = (ui)-1;
    //     goto VESREFINE_EXIT;
    // }

    // update valid_cans to index_cans, which is used for edges(index.index_)
    influenced[u] = false;
    index.index_cans_[u] = index.valid_cans_[u].top();
    index.index_cnt_[u] = index.valid_cnt_[u].top();
    cur_inf.emplace_back(u);
    for (ui i = 0; i < q_num; i++) {
        if (influenced[i]) {
            index.index_cans_[i] = index.valid_cans_[i].top();
            index.index_cnt_[i] = index.valid_cnt_[i].top();
            cur_inf.emplace_back(i);
        }
    }
    index.influenced_u_.emplace(std::move(cur_inf));

    bsxRefine_EXIT:
    delete[] influenced;
    delete[] influenced_cans;
    return returned_value;
}

void
EvaluateQuery::bsxGenResult(ui indep_num, const VertexID* indep, BSXIndex& index) {
    auto& visited_v = index.visited_v;
    auto& indep_con_cnt = index.indep_con_cnt_;
    auto& sep_flag = index.sep_flag_;  // indexed by idx
    auto& embedding_cnt = index.level_embeddings_;
    auto& label_embeddings = index.label_embeddings_;
    auto qnum = index.q_graph_->getVerticesCount();
    auto label_num = index.q_graph_->getLabelsCount();
    // generate valid_cans of indeps
    std::vector<std::vector<VertexID>> cans;
    cans.resize(qnum);
    for (ui i = 0; i < index.num_cover_; i++) {
        auto u = index.embedding->depth2u[i];
        cans[u].reserve(index.valid_cnt_[u].top());
        for (ui j = 0; j < index.valid_cnt_[u].top(); j++) {
            cans[u].emplace_back(index.valid_cans_[u].top()[j]);
        }
    }
    if (bsxGenIndepValidCans(indep_num, indep, index, cans) == false) return;
    mpz_set_ui(embedding_cnt, 1);
    for (ui l_idx = 0; l_idx < label_num; l_idx++) {
        ui nodes_num;
        // these nodes have the same label
        auto nodes = index.q_graph_->getVerticesByLabel(l_idx, nodes_num);
        if (nodes_num == 0) continue;
        // compute the number of valid embedding
        // 1.By intersected, compute the cans which may conflict with others
        // 2.Based on conflict info, seperate cans into two part
        //   con: may conflict with others nodes, process as a backtracking
        //   un-con: will not conflict with others, use (#un-con) * (#embeddings of down-level)
        //   con&un-con->all need process visited_v
        // 3.Order does not matter in this backtracking, and will not influence up-level
        // ** because there is no great idea to process up-coflict nodes,
        //    we do not seperate nodes base on up-conlict
        if (nodes_num == 1) {  // the order of indep are always the same
            mpz_mul_ui(embedding_cnt, embedding_cnt, cans[nodes[0]].size());
            continue;
        }
        // 1.first scan, compute upward conflict
        for (ui i = 0; i < nodes_num; i++) {
            auto& node = nodes[i];
            auto& v_cans = cans[node];
            // do not seperate nodes base on up-conlict
            // auto& upward_sep = sep_flag[cur_idx][1];
            // auto v_cans_cnt = v_cans.size();
            // int forward_idx = 0;
            // int backward_idx = v_cans_cnt - 1;  // backward_idx may be -1
            // upward_sep = bsxSepDiff(v_cans, indep_con_cnt, forward_idx, backward_idx);
            for (auto v_can:v_cans) indep_con_cnt[v_can]++;
        }
        // 2.second scan, compute downward conflict
        for (ui i = 0; i < nodes_num; i++) {
            auto& node = nodes[i];
            auto& v_cans = cans[node];
            auto& downward_sep0 = sep_flag[node][0];  // 0->sep the up-conflicts
            // auto& downward_sep1 = sep_flag[cur_idx][2];  // 1->sep the up-uncon.
            auto v_cans_cnt = v_cans.size();  // assert(v_cans_cnt > 0) check at cans generation
            int forward_idx = 0;
            // int middle_idx = sep_flag[cur_idx][1];
            // int backward_idx = v_cans_cnt - 1;
            for (auto v_can:v_cans) indep_con_cnt[v_can]--;
            downward_sep0 = bsxSepDiff(v_cans, indep_con_cnt, forward_idx, v_cans_cnt - 1);
            // downward_sep1 = bsxSepDiff(v_cans, indep_con_cnt, middle_idx, backward_idx);
        }
        // 3.enumerate the nodes based on diff features of 4 parts
        // ** just 2 parts so far
        bsxEnumerate4Parts(sep_flag, nodes, nodes_num, cans, visited_v, label_embeddings);
        mpz_mul(embedding_cnt, embedding_cnt, label_embeddings);
        if (mpz_cmp_ui(embedding_cnt, 0) == 0) return;
    }

    return;
}

// according to indep_con_cnt info, seperate v_cans into two parts, return the #first_part(true)
ui
EvaluateQuery::bsxSepDiff(std::vector<VertexID> &v_cans, const ui *indep_con_cnt, int forward_idx, int backward_idx) {
    if (backward_idx-forward_idx == 0) return indep_con_cnt[v_cans[forward_idx]] != 0;
    ui first_con_cnt = indep_con_cnt[v_cans[forward_idx]];
    VertexID first_idx = v_cans[forward_idx];
    while(forward_idx < backward_idx) {
        while(forward_idx < backward_idx && !indep_con_cnt[v_cans[backward_idx]]) backward_idx--;
        if (forward_idx < backward_idx)
            v_cans[forward_idx++] = v_cans[backward_idx];
        while(forward_idx < backward_idx && indep_con_cnt[v_cans[forward_idx]]) forward_idx++;
        if (forward_idx < backward_idx)
            v_cans[backward_idx--] = v_cans[forward_idx];
    }
    v_cans[forward_idx] = first_idx;
    if (first_con_cnt) forward_idx++;
    return forward_idx;
}

// TODO: opt to three parts
void  // 4 parts: up-down,up-x,x-down,x-x; down&x 2 parts so far
EvaluateQuery::bsxEnumerate4Parts(ui **&sep_flags, const VertexID* nodes, ui nodes_num,
                                   std::vector<std::vector<VertexID>>& cans, bool *&visited_v,
                                   mpz_t cur_cnt) {
    ui depth = 0;
    ui* idx = new ui[nodes_num];
    ui* cnt = new ui[nodes_num];
    ui* un_con_cnt = new ui[nodes_num];
    idx[depth] = 0;
    cnt[depth] = sep_flags[nodes[depth]][0] < cans[nodes[depth]].size()
                 ? sep_flags[nodes[depth]][0] + 1 : cans[nodes[depth]].size();
    mpz_t* embedding_level = new mpz_t[nodes_num];
    for (ui i = 0; i < nodes_num; i++) {
        mpz_init(embedding_level[i]);
    }
    mpz_set_ui(embedding_level[depth], 0);
    un_con_cnt[depth] = 0;
    while (true) {
        while (idx[depth] < cnt[depth]) {
            auto u_idx = nodes[depth];
            VertexID& v = cans[u_idx][idx[depth]];
            ui& cur_sep = sep_flags[u_idx][0];
            if (depth == nodes_num - 1) {
                ui tmp_cnt = 0;
                for (ui i = 0; i < cans[u_idx].size(); i++) {
                    if (!visited_v[cans[u_idx][i]]) tmp_cnt++;
                }
                mpz_add_ui(embedding_level[depth], embedding_level[depth], tmp_cnt);
                break;
            } else {
                idx[depth]++;
                if (idx[depth] > cur_sep) {
                    for (ui i = cur_sep; i < cans[u_idx].size(); i++) {
                        auto& can = cans[u_idx][i];
                        if (!visited_v[can]) un_con_cnt[depth]++;
                    }
                    if (un_con_cnt[depth] == 0) break;
                } else {
                    if (visited_v[v]) continue;
                    visited_v[v] = true;
                }
                depth++;
                idx[depth] = 0;
                cnt[depth] = sep_flags[nodes[depth]][0] + 1;
                mpz_set_ui(embedding_level[depth], 0);
                un_con_cnt[depth] = 0;
            }
        }
        depth--;
        if (depth == (ui)-1) {
            break;
        }
        if (idx[depth] > sep_flags[nodes[depth]][0]) {
            // process nodes which will not conflict downward
            mpz_mul_ui(embedding_level[depth+1], embedding_level[depth+1], un_con_cnt[depth]);
        } else {
            // process visited_v
            VertexID& v = cans[nodes[depth]][idx[depth]-1];
            visited_v[v] = false;
        }
        mpz_add(embedding_level[depth], embedding_level[depth], embedding_level[depth+1]);
    }

    delete[] idx;
    delete[] cnt;
    delete[] un_con_cnt;
    mpz_set(cur_cnt, embedding_level[0]);
    for (ui i = 0; i < nodes_num; i++) {
        mpz_clear(embedding_level[i]);
    }
    delete[] embedding_level;
    return;
}
