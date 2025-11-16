#include "EvaluateQuery.h"
#include "computesetintersection.h"
#include "rapidMatch/execution_tree/execution_tree_generator.h"
#include "bsx/IndepSet.h"
#include "bsx/nodeSim.h"
#include "timeOp.h"
#include <vector>
#include <cstring>

#include "pretty_print.h"

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

bool
EvaluateQuery::BS1Engine(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
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
                exploreGenValidCanIdx(data_graph, cur_depth, embedding, idx_embedding, idx_count,
                                            valid_candidate_idx, edge_matrix, visited_vertices, bn,
                                            bn_count, order, pivot, candidates, query_graph);
            }
        }

        // backtrack
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

void EvaluateQuery::exploreGenValidCanIdx(const Graph *data_graph, ui depth, ui *embedding, ui *idx_embedding,
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

void EvaluateQuery::kssComValidCans(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count,
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
    ui kernel_num = 0, shell_num = 0;
    ui* kernel = new ui[q_num];
    ui* shell = new ui[q_num];
    // kernel(true) or shell(false)
    bool* kos = new bool[q_num];
    memset(kos, 0, sizeof(bool)*q_num);
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
            degree[nbrs[j]]--;
        }
    }

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

    std::vector<VertexID> update;

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = kernel[cur_depth];
            VertexID v = valid_cans[u][idx[cur_depth]];

            embedding[u] = v;
            visited_v[v] = true;
            idx[cur_depth] += 1;

            if (TimeOp::getClockNan() >= time_limit) {
                overtime = true;
                goto EXIT;
            }

            updateShell2Kernel(query_graph, u, shell2kernel, kos, update);

            if (cur_depth == kernel_num - 1) {
                for(ui i = 0; i < shell_num; i++) {
                    VertexID u_shell = shell[i];
                    kssComValidCans(data_graph, query_graph, candidates, candidates_count, valid_cans,
                                 valid_cans_cnt, embedding, u_shell, visited_u, visited_v);
                }
                if (kssGenResult(shell_num, shell, valid_cans, valid_cans_cnt, visited_v, embedding_cnt, time_limit) == true) {
                    overtime = true;
                    goto EXIT;
                }
                visited_v[v] = false;
                if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
                    goto EXIT;
                }
                restoreShell2Kernel(query_graph, u, shell2kernel, kos);
            } else {
                cur_depth += 1;
                VertexID next_u = kernel[cur_depth];
                idx[cur_depth] = 0;
                call_cnt += 1;
                kssComValidCans(data_graph, query_graph, candidates, candidates_count, valid_cans,
                                 valid_cans_cnt, embedding, next_u, visited_u, visited_v);
                
                visited_u[next_u] = true;
                idx_count[cur_depth] = valid_cans_cnt[next_u];
            }
        }

        cur_depth -= 1;
        if (cur_depth == (ui)-1)
            break;
        VertexID last_u = kernel[cur_depth + 1];
        visited_v[embedding[kernel[cur_depth]]] = false;
        visited_u[last_u] = false;
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

void
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

void
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

bool 
EvaluateQuery::kssGenResult(ui shell_num, ui* shell, ui** valid_cans, ui* valid_cans_count, bool * visited_v,
                                        mpz_t embedding_cnt, int64_t& time_limit) {
    
    return kssGenResultImpl(0, shell_num, shell, valid_cans, valid_cans_count, visited_v, embedding_cnt, time_limit);
}

bool 
EvaluateQuery::kssGenResultImpl(ui depth, ui shell_num, ui* shell, ui** valid_cans, ui* valid_cans_count, bool * visited_v,
                                            mpz_t embedding_cnt, int64_t& time_limit) {

    if (TimeOp::getClockNan() >= time_limit) {
        return true;
    }
    VertexID u_shell = shell[depth];

    if (depth == shell_num - 1) {
        for(ui i = 0; i < valid_cans_count[u_shell]; i++) {
            VertexID v_id = valid_cans[u_shell][i];
            if (!visited_v[v_id]) mpz_add_ui(embedding_cnt, embedding_cnt, 1);
        }
    }
    else {
        for(ui i = 0; i < valid_cans_count[u_shell]; i++) {
            VertexID v_id = valid_cans[u_shell][i];
            if(!visited_v[v_id]){
                visited_v[v_id] = true;
                if (kssGenResultImpl(depth+1, shell_num, shell, valid_cans, valid_cans_count, visited_v,
                                                 embedding_cnt, time_limit) == true)
                    return true;
                visited_v[v_id] = false;
            }
        }
    }
    return false;
}

// int64_t BSXIndex::getNeighbors_time = 0;
// int64_t BSXIndex::update_time = 0;
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
    // int64_t enumerate_time = 0;
    // int64_t batch_time = 0;
    // int64_t refine_time = 0;
    // int64_t start_time;
    ui* order = nullptr;
    // separate leaf and trunk vertices(min_vertex_cover)
    ui num_cover = 0;
    maxCoverOrder(query_graph, order, num_cover, candidates_count);
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

    // init info of start vertex
    batch_info[start_vertex].add();
    // start_time = TimeOp::getClockNan();
    bsxComEqBatch(index, start_vertex);
    // batch_time += TimeOp::getClockNan() - start_time;
    index.valid_cans_[start_vertex].push(new VertexID[batch_info[start_vertex].maxCnt_.top()]);
    index.valid_cnt_[start_vertex].push(0);

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

    if (status == 0) {
        std::cout << g_name << ": Folder created successfully.\n";
    } else {
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

            // start_time = TimeOp::getClockNan();
            VertexID failed_u = bsxRefine(index, u);
            // refine_time += TimeOp::getClockNan() - start_time;
            if (failed_u != (ui)-1) {  // no valid cans for next depth
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
                // start_time = TimeOp::getClockNan();
                bsxGenResult(num_indep, indep_nodes, index);
                // enumerate_time += TimeOp::getClockNan() - start_time;
                mpz_add(embedding_cnt, embedding_cnt, level_embeddings);
                if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
                    goto EXIT;
                }
                // next batch
                bsxDeRefine(index);
            } else {
                cur_depth++;
                VertexID cur_u = bsxGenNxtU(index, order, cur_depth, num_cover);
                if (cur_u == (VertexID)-1) cur_u = order[cur_depth];
                depth2u.emplace_back(cur_u);
                call_count++;
                // construct nbrs&seperate batches, and then refinement
                batch_info[cur_u].add();
                // start_time = TimeOp::getClockNan();
                bsxComEqBatch(index, cur_u);
                // batch_time += TimeOp::getClockNan() - start_time;
                index.valid_cans_[cur_u].push(new VertexID[batch_info[cur_u].maxCnt_.top()]);
                index.valid_cnt_[cur_u].push(0);
                visited_u[cur_u] = true;
            }
        }

        // backtracking
        cur_depth -= 1;
        if (cur_depth == ui(-1))
            break;
        VertexID last_u = depth2u[cur_depth+1];
        depth2u.resize(cur_depth+1);
        visited_u[last_u] = false;

        batch_info[last_u].pop();
        delete[] index.valid_cans_[last_u].top();
        index.valid_cans_[last_u].pop();
        index.valid_cnt_[last_u].pop();
        bsxDeRefine(index);
    }

    // Release the buffer.
    EXIT:

#ifdef ANALYZE_DUPLICATE
    for (ui i = 1; i < q_num; i++) {
        out_files[i-1].close();
    }
#endif
    // std::cout << "enumerate_time: " << enumerate_time << std::endl;
    // std::cout << "refine_time: " << refine_time << std::endl;
    // std::cout << "batch_time: " << batch_time << std::endl;
    // std::cout << "getNeighbors_time: " << BSXIndex::getNeighbors_time << std::endl;
    // std::cout << "update_time: " << BSXIndex::update_time << std::endl;
    return overtime;
}

/**
 * generate order for BSXEngine
 * 1. static: seperate vertices into cover&independent by MaxCover
 * 2. dynamic: sort each kind by #cans(asc), #degree(des), #id(asc)
 * compute static order by max indep cover, dynamic order computed along matching
*/
void
EvaluateQuery::maxCoverOrder(const Graph *graph, ui *&order, ui& num_cover, ui *candidates_count) {
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
    // memset(nbr_updated, false, sizeof(bool)*q_num);
    std::copy(index.visited_u, index.visited_u+q_num, nbr_updated);

    // process first u in influenced_u seperately, valid_cans of first_influenced_u comes from
    //   its batch_nodes, and couldn't be deleted
    auto first_u = index.influenced_u_.top()[0];
    // nbr_updated[first_u] = true;
    // remove index_, added at refinement
    ui nbrs_cnt;
    auto nbrs = index.q_graph_->getVertexNeighbors(first_u, nbrs_cnt);
    for (ui i = 0; i < nbrs_cnt; i++) {
        auto& nbr = nbrs[i];
        if (nbr_updated[nbr]) continue;
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

// equ-batch refine, just process first v of valid_cans, because of they are equ
// if success, return -1, else return failed uId
ui
EvaluateQuery::bsxRefine(BSXIndex& index, VertexID u) {
    // int64_t start_time;
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
    // start_time = TimeOp::getClockNan();
    index.updateIndex(influenced, u);
    // BSXIndex::update_time += TimeOp::getClockNan() - start_time;
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
            // upward_sep = sepDiff(v_cans, indep_con_cnt, forward_idx, backward_idx);
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
            downward_sep0 = sepDiff(v_cans, indep_con_cnt, forward_idx, v_cans_cnt - 1);
            // downward_sep1 = sepDiff(v_cans, indep_con_cnt, middle_idx, backward_idx);
        }
        // 3.enumerate the nodes based on diff features of 4 parts
        // ** just 2 parts so far
        enum4Parts(sep_flag, nodes, nodes_num, cans, visited_v, label_embeddings);
        mpz_mul(embedding_cnt, embedding_cnt, label_embeddings);
        if (mpz_cmp_ui(embedding_cnt, 0) == 0) return;
    }

    return;
}

// according to indep_con_cnt info, seperate v_cans into two parts, return the #first_part(true)
ui
EvaluateQuery::sepDiff(std::vector<VertexID> &v_cans, const ui *indep_con_cnt, int forward_idx, int backward_idx) {
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
EvaluateQuery::enum4Parts(ui **&sep_flags, const VertexID* nodes, ui nodes_num,
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

vector<ui> EdgeSub::down_record;
/**
 * use FiPE method
*/
bool
EvaluateQuery::FiPEEngine(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix,
                             ui **candidates, ui *candidates_count,
                             size_t output_limit_num, size_t &call_count, mpz_t embedding_cnt, int64_t& time_limit) {
    auto qnum = query_graph->getVerticesCount();
    VertexID* order = new VertexID[qnum];
    ui num_cover = FiPEIndep::indepSetOnDegree(query_graph, order);
    // main data structure
    FiPEIndex index(query_graph, data_graph, edge_matrix, candidates, candidates_count, num_cover, order);
    EdgeSub::down_record.resize(data_graph->getVerticesCount(), (ui)-1);

    // auxiliary data structure
    bool overtime = false;
    mpz_init_set_ui(embedding_cnt, 0);
    ui cur_depth = 0;
    VertexID start_vertex = order[cur_depth];
    auto& level_embeddings = index.indepInfo->embedding_total;
    auto& subInfo = index.subInfo_;

    splitCans(index, cur_depth);
    while (!comStartCans(index) || !comCurSpace(index, cur_depth)) {
        if (nxtSubCans(index, cur_depth) == false) goto EXIT;
    }
    comSub(index, cur_depth);

    while (true) {
        while (subInfo[cur_depth].subs.cur_s < subInfo[cur_depth].subs.c_cnt) {
            VertexID u = order[cur_depth];

            if (TimeOp::getClockNan() >= time_limit) {
                overtime = true;
                goto EXIT;
            }

            if (subInfo[cur_depth].connected) setCurSpaceCon(index, cur_depth);
            else setCurSpaceDis(index, cur_depth);
            subInfo[cur_depth].subs.cur_s++;

            if (cur_depth >= num_cover - 2) {
                // enumerate results on indep nodes, process ancestors' ves by the way
                // mpz_ui_sub(index.indepInfo->remained, output_limit_num, embedding_cnt);
                FiPEEnum(index);
                mpz_add(embedding_cnt, embedding_cnt, level_embeddings);
                if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
                    goto EXIT;
                }
                // nxt batch
                if (subInfo[cur_depth].connected) clearCurSpaceCon(index, cur_depth);
                else clearCurSpaceDis(index, cur_depth);
            } else {
                cur_depth++;
                splitCans(index, cur_depth);
                bool nxt_valid = true;
                while (!comCurSpace(index, cur_depth)) {
                    if (!nxtSubCans(index, cur_depth)) {
                        index.subCans_[cur_depth+1].splitted = false;
                        nxt_valid = false;
                        break;
                    }
                }
                if (!nxt_valid) break;
                comSub(index, cur_depth);
                call_count++;
            }
        }

        // backtracking
        bool nxt_valid = true;
        if (nxtSubCans(index, cur_depth)) {
            if (cur_depth == 0) {
                while (!comStartCans(index) || !comCurSpace(index, cur_depth)) {
                    if (!nxtSubCans(index, cur_depth)) {
                        nxt_valid = false;
                        break;
                    }
                }
            } else {
               while (!comCurSpace(index, cur_depth)) {
                    if (nxtSubCans(index, cur_depth) == false) {
                        nxt_valid = false;
                        break;
                    }
                }
            }
            if (nxt_valid) {
                comSub(index, cur_depth);
                continue;
            }
        }
        cur_depth--;
        if (cur_depth == ui(-1))
            break;
        if (subInfo[cur_depth].connected) clearCurSpaceCon(index, cur_depth);
        else clearCurSpaceDis(index, cur_depth);
    }

    // Release the buffer.
    EXIT:
    return overtime;
}


/**
 * split cans when there are too much cans
 */
inline void
EvaluateQuery::splitCans(FiPEIndex& index, ui depth) {
    auto& down = index.order_[depth+1];
    auto& subCans = index.subCans_[depth+1];
    subCans.up_changed = true;
    if (depth == 0) {
        auto& up = index.order_[depth];
        auto& upCans = index.subCans_[depth];
        if (index.valid_cans_[up].back().size() > MIN_SUBCANS) {
            upCans.splitted = true;
            upCans.cur_start = 0;
            upCans.cur_end = MIN_SUBCANS;
            vector<VertexID> new_cans;
            new_cans.reserve(MIN_SUBCANS);
            new_cans.insert(new_cans.end(),
                            index.valid_cans_[up].back().begin(),
                            index.valid_cans_[up].back().begin()+MIN_SUBCANS);
            index.valid_cans_[up].push_back(move(new_cans));
        } else {
            upCans.splitted = false;
        }
    }
    // if too much candidates, split
    if (!index.subInfo_[depth].connected
        && index.valid_cans_[down].back().size() > MIN_SUBCANS) {
        subCans.splitted = true;
        subCans.cur_start = 0;
        subCans.cur_end = MIN_SUBCANS;
        vector<VertexID> new_cans;
        new_cans.reserve(MIN_SUBCANS);
        new_cans.insert(new_cans.end(),
                        index.valid_cans_[down].back().begin(),
                        index.valid_cans_[down].back().begin()+MIN_SUBCANS);
        index.valid_cans_[down].push_back(move(new_cans));
    } else {
        subCans.splitted = false;
    }
    return;
}

/**
 * nxt sub cans
 * false->splited but have no next subset of cans, true->o.w.
 */
inline bool
EvaluateQuery::nxtSubCans(FiPEIndex& index, ui depth) {
    auto& down = index.order_[depth+1];
    auto& subCans = index.subCans_[depth+1];
    subCans.up_changed = false;
    if (subCans.splitted) index.valid_cans_[down].pop_back();
    if (depth == 0) {
        auto& up = index.order_[depth];
        auto& upCans = index.subCans_[depth];
        if (subCans.splitted == false
            || subCans.cur_end == index.valid_cans_[down].back().size()) {
            if (upCans.splitted) index.valid_cans_[up].pop_back();
            if (upCans.splitted == false
                || upCans.cur_end == index.valid_cans_[up].back().size()) {
                return false;
            }
            upCans.cur_start = upCans.cur_end;
            upCans.cur_end = min(upCans.cur_end + MIN_SUBCANS,
                                 (ui)index.valid_cans_[up].back().size());
            vector<VertexID> new_cans;
            new_cans.reserve(MIN_SUBCANS);
            new_cans.insert(new_cans.end(),
                            index.valid_cans_[up].back().begin()+upCans.cur_start,
                            index.valid_cans_[up].back().begin()+upCans.cur_end);
            index.valid_cans_[up].push_back(move(new_cans));
            subCans.up_changed = true;
            subCans.cur_start = 0;
            subCans.cur_end = MIN_SUBCANS;
        } else {
            subCans.cur_start = subCans.cur_end;
            subCans.cur_end = min(subCans.cur_end + MIN_SUBCANS,
                                  (ui)index.valid_cans_[down].back().size());
        }
    } else {
        if (subCans.splitted == false
            || subCans.cur_end == index.valid_cans_[down].back().size())
            return false;
        subCans.cur_start = subCans.cur_end;
        subCans.cur_end = min(subCans.cur_end + MIN_SUBCANS,
                              (ui)index.valid_cans_[down].back().size());
    }
    if (subCans.splitted) {
        vector<VertexID> new_cans;
        new_cans.reserve(MIN_SUBCANS);
        new_cans.insert(new_cans.end(),
                        index.valid_cans_[down].back().begin()+subCans.cur_start,
                        index.valid_cans_[down].back().begin()+subCans.cur_end);
        index.valid_cans_[down].push_back(move(new_cans));
    }
    return true;
}

/**
 * 1.set valid_cans of cur_s
 * 3.push valid_cans from inf_cans to FiPEIndex
 * 1.generate sub_cans for next_u
 *   1.1 connected: union of nbrs from (sub_cans of u) to next_u
 *   1.2 o.w.: u will not influence the cans of next_u
 * 2.compute the valid edges for current u & u_n-1(if conncected)
 * return: false, if valid_cans of delayed_nbrs is empty
*/
bool
EvaluateQuery::setCurSpaceCon(FiPEIndex& index, ui depth) {
    auto& subInfo = index.subInfo_[depth];
    auto& subs = subInfo.subs;
    auto& down_cans = subInfo.down_cans;  // down_idx->down_cans
    auto& down_idxs = subInfo.down_idxs;  // edge_idx->down_idx
    auto& edge_up_idxs = subs.edge_up_idxs;
    auto& edge_down_idxs = subs.edge_down_idxs;
    auto& cur_s = subs.cur_s;
    auto& edge_up_start = subs.edge_up_offset[cur_s];
    auto& edge_up_end = subs.edge_up_offset[cur_s+1];
    auto& edge_down_start = subs.edge_down_offset[edge_up_start];
    auto& edge_down_end = subs.edge_down_offset[edge_up_end];
    auto& down_record = subs.down_record;
    auto& influenced = subInfo.influenced;
    auto& inf_cans = subInfo.inf_cans;
    auto& valid_cans = index.valid_cans_;

    if (depth < index.indepInfo->num_cover_-2) {
        // compute up_cans of nxtSub
        auto& nxtSub = index.subInfo_[depth+1];
        auto& nxt_up_cans = nxtSub.up_cans;
        nxt_up_cans.clear();
        auto& nxt_up_idxs = subInfo.nxt_up_idxs;
        auto& last_down_idxs = nxtSub.last_down_idxs;
        nxt_up_idxs.resize(down_cans.size());
        last_down_idxs.clear();
        for (ui i = edge_down_start; i < edge_down_end; i++) {
            auto& edge_idx = edge_down_idxs[i];
            auto& down_idx = down_idxs[edge_idx];  // idx of influenced
            auto& down = down_cans[down_idx];
            if (down_record[down] == (ui)-1) {
                down_record[down] = nxt_up_cans.size();
                nxt_up_cans.emplace_back(down);
                last_down_idxs.emplace_back(down_idx);
                nxt_up_idxs[down_idx] = down_record[down];
            }
        }
        for (auto& nxt_up_can : nxt_up_cans) down_record[nxt_up_can] = (ui)-1;

        // compute the valid_cans of delayed_nbrs
        auto& delayed_nbrs = subInfo.nbrs.delayed_;
        subInfo.delayed_inf.clear();
        subInfo.delayed_inf.resize(delayed_nbrs.size(), true);
        for (ui i = 0; i < delayed_nbrs.size(); i++) {
            auto& delayed_nbr = delayed_nbrs[i];
            // if any can do not influence shared_nbr, valid_cans keeps
            for (ui j = 0; j < nxt_up_cans.size(); j++) {
                auto& down_idx = last_down_idxs[j];
                if (!influenced[delayed_nbr][down_idx]) {
                    subInfo.delayed_inf[i] = false;
                    break;
                }
            }
            if (subInfo.delayed_inf[i]) {
                vector<vector<VertexID>> all_cans;
                all_cans.reserve(nxt_up_cans.size());
                for (ui j = 0; j < nxt_up_cans.size(); j++) {
                    auto& down_idx = last_down_idxs[j];
                    all_cans.emplace_back(inf_cans[delayed_nbr][down_idx]);
                }
                valid_cans[delayed_nbr].push_back(move(SetOp::unionMultiple(all_cans)));
            }
        }
    }

    // push valid_cans from inf_cans to FiPEIndex
    // use the first as the representative
    auto& up_idx = edge_up_idxs[edge_up_start];
    for (auto& nbr : subInfo.nbrs.up_indep_) {
        if (influenced[nbr][up_idx]) {
            valid_cans[nbr].push_back(inf_cans[nbr][up_idx]);
        }
    }
    auto& edge_idx = edge_down_idxs[edge_down_start];
    for (auto& nbr : subInfo.nbrs.shared_) {
        if (influenced[nbr][edge_idx]) {
            valid_cans[nbr].push_back(move(inf_cans[nbr][edge_idx]));
        }
    }
    auto down_idx = down_idxs[edge_idx];
    for (auto& nbr : subInfo.nbrs.down_indep_) {
        if (influenced[nbr][down_idx]) {
            valid_cans[nbr].push_back(inf_cans[nbr][down_idx]);
        }
    }
    return true;
}

/**pop the valid_cans of influenced_u from FiPEIndex
*/
void
EvaluateQuery::clearCurSpaceCon(FiPEIndex& index, ui depth) {
    auto& subInfo = index.subInfo_[depth];
    auto& subs = subInfo.subs;
    auto& down_idxs = subInfo.down_idxs;  // edge_idx->down_idx
    auto& edge_up_idxs = subs.edge_up_idxs;
    auto& edge_down_idxs = subs.edge_down_idxs;
    auto& cur_s = subs.cur_s;
    auto& edge_up_start = subs.edge_up_offset[cur_s-1];
    auto& edge_down_start = subs.edge_down_offset[edge_up_start];
    auto& influenced = subInfo.influenced;
    auto& valid_cans = index.valid_cans_;
    auto& delayed_inf = subInfo.delayed_inf;
    // pop valid_cans from inf_cans to FiPEIndex
    // use the first as the representative
    auto& up_idx = edge_up_idxs[edge_up_start];
    for (auto& nbr : subInfo.nbrs.up_indep_) {
        if (influenced[nbr][up_idx]) {
            valid_cans[nbr].pop_back();
        }
    }
    auto& edge_idx = edge_down_idxs[edge_down_start];
    for (auto& nbr : subInfo.nbrs.shared_) {
        if (influenced[nbr][edge_idx]) {
            valid_cans[nbr].pop_back();
        }
    }
    auto down_idx = down_idxs[edge_idx];
    for (auto& nbr : subInfo.nbrs.down_indep_) {
        if (influenced[nbr][down_idx]) {
            valid_cans[nbr].pop_back();
        }
    }
    for (ui i = 0; i < subInfo.nbrs.delayed_.size(); i++) {
        if (delayed_inf[i]) {
            valid_cans[subInfo.nbrs.delayed_[i]].pop_back();
        }
    }
}

void
EvaluateQuery::setCurSpaceDis(FiPEIndex& index, ui depth) {
    auto& subInfo = index.subInfo_[depth];
    auto& subs = subInfo.subs;
    auto& down_cans = subInfo.down_cans;  // down_idx->down_cans
    auto& edge_down_idxs = subs.edge_down_idxs;
    auto& cur_s = subs.cur_s;
    auto& edge_down_start = subs.edge_down_offset[cur_s];
    auto& edge_down_end = subs.edge_down_offset[cur_s+1];
    auto& down_record = subs.down_record;
    auto& influenced = subInfo.influenced;
    auto& inf_cans = subInfo.inf_cans;
    auto& valid_cans = index.valid_cans_;
    ui down_num = down_cans.size();

    if (depth < index.indepInfo->num_cover_-2) {
        // compute up_cans of nxtSub
        auto& nxtSub = index.subInfo_[depth+1];
        auto& nxt_up_cans = nxtSub.up_cans;
        auto& nxt_up_idxs = subInfo.nxt_up_idxs;
        auto& last_down_idxs = nxtSub.last_down_idxs;
        nxt_up_idxs.resize(down_cans.size());
        last_down_idxs.clear();
        nxt_up_cans.clear();
        for (ui i = edge_down_start; i < edge_down_end; i++) {
            auto& edge_idx = edge_down_idxs[i];  // idx of influenced(shared)
            ui down_idx = edge_idx%down_num;       // idx of influenced(indep)
            auto& down = down_cans[down_idx];
            if (down_record[down] == (ui)-1) {
                down_record[down] = nxt_up_cans.size();
                nxt_up_cans.emplace_back(down);
                last_down_idxs.emplace_back(down_idx);
                nxt_up_idxs[down_idx] = down_record[down];
            }
        }
        for (auto& nxt_up_can : nxt_up_cans) down_record[nxt_up_can] = (ui)-1;

        // compute the valid_cans of delayed_nbrs
        auto& delayed_nbrs = subInfo.nbrs.delayed_;
        subInfo.delayed_inf.clear();
        subInfo.delayed_inf.resize(delayed_nbrs.size(), true);
        for (ui i = 0; i < delayed_nbrs.size(); i++) {
            auto& delayed_nbr = delayed_nbrs[i];
            // if any can do not influence shared_nbr, valid_cans keeps
            subInfo.delayed_inf[i] = true;
            for (ui j = 0; j < nxt_up_cans.size(); j++) {
                auto& down_idx = last_down_idxs[j];
                if (!influenced[delayed_nbr][down_idx]) {
                    subInfo.delayed_inf[i] = false;
                    break;
                }
            }
            if (subInfo.delayed_inf[i]) {
                vector<vector<VertexID>> all_cans;
                all_cans.reserve(nxt_up_cans.size());
                for (ui j = 0; j < nxt_up_cans.size(); j++) {
                    auto& down_idx = last_down_idxs[j];
                    all_cans.emplace_back(inf_cans[delayed_nbr][down_idx]);
                }
                valid_cans[delayed_nbr].push_back(move(SetOp::unionMultiple(all_cans)));
            }
        }
    }

    // push valid_cans from inf_cans to FiPEIndex
    // use the first as the representative
    auto& edge_idx = edge_down_idxs[edge_down_start];
    auto up_idx = edge_idx/down_num;
    auto down_idx = edge_idx%down_num;
    for (auto& nbr : subInfo.nbrs.up_indep_) {
        if (influenced[nbr][up_idx]) {
            valid_cans[nbr].push_back(inf_cans[nbr][up_idx]);
        }
    }
    for (auto& nbr : subInfo.nbrs.shared_) {
        if (influenced[nbr][edge_idx]) {
            valid_cans[nbr].push_back(move(inf_cans[nbr][edge_idx]));
        }
    }
    for (auto& nbr : subInfo.nbrs.down_indep_) {
        if (influenced[nbr][down_idx]) {
            valid_cans[nbr].push_back(inf_cans[nbr][down_idx]);
        }
    }
}

/**pop the valid_cans of influenced_u from FiPEIndex
*/
void
EvaluateQuery::clearCurSpaceDis(FiPEIndex& index, ui depth) {
    auto& subInfo = index.subInfo_[depth];
    auto& subs = subInfo.subs;
    auto& down_cans = subInfo.down_cans;
    auto& edge_down_idxs = subs.edge_down_idxs;
    auto& cur_s = subs.cur_s;
    auto& edge_down_start = subs.edge_down_offset[cur_s-1];
    auto& influenced = subInfo.influenced;
    auto& valid_cans = index.valid_cans_;
    auto& delayed_inf = subInfo.delayed_inf;
    ui down_num = down_cans.size();
    // use the first as the representative
    auto& edge_idx = edge_down_idxs[edge_down_start];
    auto up_idx = edge_idx/down_num;
    auto down_idx = edge_idx%down_num;
    for (auto& nbr : subInfo.nbrs.up_indep_) {
        if (influenced[nbr][up_idx]) {
            valid_cans[nbr].pop_back();
        }
    }
    for (auto& nbr : subInfo.nbrs.shared_) {
        if (influenced[nbr][edge_idx]) {
            valid_cans[nbr].pop_back();
        }
    }
    for (auto& nbr : subInfo.nbrs.down_indep_) {
        if (influenced[nbr][down_idx]) {
            valid_cans[nbr].pop_back();
        }
    }
    for (ui i = 0; i < subInfo.nbrs.delayed_.size(); i++) {
        if (delayed_inf[i]) {
            valid_cans[subInfo.nbrs.delayed_[i]].pop_back();
        }
    }
}

void
EvaluateQuery::FiPEEnum(FiPEIndex& index) {
#ifdef ANALYZE_TIME
    auto start_time = TimeOp::getClockNan();
#endif
    auto qnum = index.q_graph_->getVerticesCount();
    auto dnum = index.d_graph_->getVerticesCount();
    auto& order = index.order_;
    auto& used_cans = index.indepInfo->used_cans_;
    memset(used_cans, false, sizeof(bool)*dnum);
    auto& embedding_total = index.indepInfo->embedding_total;
    mpz_set_ui(embedding_total, 0);
    auto& embedding_uncon = index.indepInfo->embedding_uncon;
    auto& uncon = index.indepInfo->uncon;
    uncon = false;
    auto& embedding_step = index.indepInfo->embedding_step;
    auto& cans = index.indepInfo->cans;
    auto& visited_v = index.visited_v;
    auto& subInfo = index.subInfo_;

    // extract the candidates of indep vertices
    for (VertexID i = 0; i < qnum; i++) {
        if (index.indepInfo->indep_bool[i]) {
            cans[i] = index.valid_cans_[i].back();
            for (auto& can:cans[i]) used_cans[can] = true;
        }
    }

    auto& num_cover = index.indepInfo->num_cover_;
    static vector<ui> idx(num_cover);
    static vector<pair<ui, ui>> up_down_idxs(num_cover-1);  // up_idx & down_idx of matched edges
    static ui depth;
    static vector<vector<ui>> edges_idxs(num_cover-1);
    depth = 0;
    ui edge_start, edge_end;
    if (subInfo[depth].connected) {
        auto& edge_up_start = subInfo[depth].subs.edge_up_offset[index.subInfo_[depth].subs.cur_s-1];
        auto& edge_up_end = subInfo[depth].subs.edge_up_offset[index.subInfo_[depth].subs.cur_s];
        edge_start = subInfo[depth].subs.edge_down_offset[edge_up_start];
        edge_end = subInfo[depth].subs.edge_down_offset[edge_up_end];
    } else {
        edge_start = subInfo[depth].subs.edge_down_offset[index.subInfo_[depth].subs.cur_s-1];
        edge_end = subInfo[depth].subs.edge_down_offset[index.subInfo_[depth].subs.cur_s];
    }
    edges_idxs[depth].clear();
    for (ui i = edge_start; i < edge_end; i++) {
        edges_idxs[depth].emplace_back(index.subInfo_[depth].subs.edge_down_idxs[i]);
    }
    idx[depth] = 0;
    static ui con_num;  // record the conflicts number
    con_num = 0;
    while(true) {
        while(idx[depth] < edges_idxs[depth].size()) {
            ui up_idx, down_idx;
            auto& edge_idx =  edges_idxs[depth][idx[depth]];
            idx[depth]++;
            if (index.subInfo_[depth].connected) {
                up_idx = index.subInfo_[depth].up_idxs[edge_idx];
                down_idx = index.subInfo_[depth].down_idxs[edge_idx];
            } else {
                auto down_num = index.subInfo_[depth].down_cans.size();
                up_idx = edge_idx/down_num;
                down_idx = edge_idx%down_num;
            }
            up_down_idxs[depth] = make_pair(up_idx, down_idx);
            auto& up = index.subInfo_[depth].up_cans[up_idx];
            auto& down = index.subInfo_[depth].down_cans[down_idx];
            if (depth == 0) {
                if (up == down) continue;
                cans[order[depth]][0] = up;
                if (used_cans[up]) con_num++;
                visited_v[up] = true;
            }
            if (visited_v[down]) continue;
            visited_v[down] = true;
            if (used_cans[down]) con_num++;
            cans[order[depth+1]][0] = down;
            if (depth >= num_cover - 2) {
                if (con_num == 0 && uncon == true) {  // no conflicts & searched the no-conflicts case
                    mpz_add(embedding_total, embedding_total, embedding_uncon);
                } else {
#ifdef FIPE_HOMOMORPHISM
                    index.indepInfo->homoEnum();
#else
                    index.indepInfo->enumeration(index.q_graph_);
#endif
                    if (con_num == 0) {  // no conflicts, record the results
                        mpz_set(embedding_uncon, embedding_step);
                        uncon = true;
                    }
                    mpz_add(embedding_total, embedding_total, embedding_step);
                }
                // if (mpz_cmp(embedding_total, index.indepInfo->remained) > 0) return;
                if (depth == 0) {
                    if (used_cans[up]) con_num--;
                    visited_v[up] = false;
                }
                if (used_cans[down]) con_num--;
                visited_v[down] = false;
            } else {
                depth++;
                // compute the valid_edges of nxt depth, based on down_idx
                if (subInfo[depth].connected) {
                    ui nxt_idx = subInfo[depth-1].nxt_up_idxs[down_idx];
                    auto& edge_up_start = subInfo[depth].subs.edge_up_offset[index.subInfo_[depth].subs.cur_s-1];
                    auto& edge_up_end = subInfo[depth].subs.edge_up_offset[index.subInfo_[depth].subs.cur_s];
                    ui edge_up_idx;
                    for (edge_up_idx = edge_up_start; edge_up_idx < edge_up_end; edge_up_idx++) {
                        if (subInfo[depth].subs.edge_up_idxs[edge_up_idx] == nxt_idx) break;
                    }
                    if (edge_up_idx == edge_up_end) {  // if nxt_idx is invalid in nxtSub
                        break;
                    }
                    edge_start = subInfo[depth].subs.edge_down_offset[edge_up_idx];
                    edge_end = subInfo[depth].subs.edge_down_offset[edge_up_idx+1];
                    edges_idxs[depth].clear();
                    for (ui i = edge_start; i < edge_end; i++) {
                        edges_idxs[depth].emplace_back(index.subInfo_[depth].subs.edge_down_idxs[i]);
                    }
                } else {
                    ui nxt_idx = subInfo[depth-1].nxt_up_idxs[down_idx];
                    auto down_num = index.subInfo_[depth].down_cans.size();
                    edge_start = subInfo[depth].subs.edge_down_offset[index.subInfo_[depth].subs.cur_s-1];
                    edge_end = subInfo[depth].subs.edge_down_offset[index.subInfo_[depth].subs.cur_s];
                    edges_idxs[depth].clear();
                    for (ui i = edge_start; i < edge_end; i++) {
                        auto& edge_idx = subInfo[depth].subs.edge_down_idxs[i];
                        if (edge_idx/down_num == nxt_idx)
                            edges_idxs[depth].emplace_back(index.subInfo_[depth].subs.edge_down_idxs[i]);
                    }
                }
                idx[depth] = 0;
            }
        }
        depth--;
        if (depth == (ui)-1) break;
        ui up_idx, down_idx;
        auto& edge_idx =  edges_idxs[depth][idx[depth]-1];
        if (index.subInfo_[depth].connected) {
            up_idx = index.subInfo_[depth].up_idxs[edge_idx];
            down_idx = index.subInfo_[depth].down_idxs[edge_idx];
        } else {
            auto down_num = index.subInfo_[depth].down_cans.size();
            up_idx = edge_idx/down_num;
            down_idx = edge_idx%down_num;
        }
        up_down_idxs[depth] = make_pair(up_idx, down_idx);
        auto& up = index.subInfo_[depth].up_cans[up_idx];
        auto& down = index.subInfo_[depth].down_cans[down_idx];
        if (depth == 0) {
            if (used_cans[up]) con_num--;
            visited_v[up] = false;
        }
        if (used_cans[down]) con_num--;
        visited_v[down] = false;
    }

#ifdef ANALYZE_TIME
    FiPEIndex::enumerate_time += TimeOp::getClockNan() - start_time;
#endif
    return;
}

/**
 * compute substitutable information for up & down cans
 */
void
EvaluateQuery::comSub(FiPEIndex& index, ui depth) {
    auto& subInfo = index.subInfo_[depth];
    auto& up_cans = subInfo.up_cans;
    auto& subs = subInfo.subs;
    subs.cur_s = 0;
    auto& up_indep = subInfo.nbrs.up_indep_;
    auto& down_indep = subInfo.nbrs.down_indep_;
    auto& influenced = subInfo.influenced;
    auto& inf_cans = subInfo.inf_cans;

    // compute sub by using 2 offset alternately
    // c_* store the complete sub info, p_* store FiPE sub info
    auto& c_cnt = subs.c_cnt;
    auto& c_offset = subs.edge_up_offset;
    auto& c_idxs = subs.edge_up_idxs;
    auto& p_cnt = subs.p_cnt;
    auto& p_offset = subs.p_edge_up_offset;
    auto& p_idxs = subs.p_edge_up_idxs;
    bool same = true;

    // init up subs
    auto& up_group_num = subInfo.subs.up_group_num;
    c_offset.resize(up_cans.size()+1);
    p_offset.resize(up_cans.size()+1);
    c_idxs.resize(up_cans.size());
    p_idxs.resize(up_cans.size());
    up_group_num.clear();
    up_group_num.resize(up_cans.size(), 0);  // 0->no group
    c_cnt = 1;
    c_offset[0] = 0; c_offset[1] = up_cans.size();
    for (ui i = 0; i < up_cans.size(); i++) c_idxs[i] = i;

    // up_indep, just separate up_cans(idxs), only work with start vertex
    for (auto& unbr : up_indep) {
        vector<bool> up_grouped(c_idxs.size(), false);
        p_cnt = 0;
        for (ui i = 0; i < c_cnt; i++) {
            auto start = c_offset[i], end = c_offset[i+1];
            for (auto idx1 = start; idx1 < end; idx1++) {
                auto& cans_idx1 = c_idxs[idx1];
                if (up_grouped[cans_idx1]) continue;
                p_cnt++;
                up_grouped[cans_idx1] = true;
                up_group_num[cans_idx1] = p_cnt;
                p_offset[p_cnt] = p_offset[p_cnt-1];
                p_idxs[p_offset[p_cnt]++] = cans_idx1;
                auto& cans1 = inf_cans[unbr][cans_idx1];
                for (auto idx2 = idx1 + 1; idx2 < end; idx2++) {
                    auto& cans_idx2 = c_idxs[idx2];
                    auto& cans2 = inf_cans[unbr][cans_idx2];
                    if (up_grouped[cans_idx2]) continue;
                    if (cans1.size() != cans2.size()) continue;
                    same = true;
                    for (ui com_idx = 0; com_idx < cans1.size(); com_idx++) {
                        if (cans1[com_idx] != cans2[com_idx]) {
                            same = false;
                            break;
                        }
                    }
                    if (same) {
                        p_idxs[p_offset[p_cnt]++] = cans_idx2;
                        up_grouped[cans_idx2] = true;
                        up_group_num[cans_idx2] = p_cnt;
                    }
                }
            }
        }
        c_idxs.swap(p_idxs);
        c_offset.swap(p_offset);
        c_cnt = p_cnt;
    }

    // down_indep, seperate down_cans(idxs)
    auto& down_group_num = subInfo.subs.down_group_num;
    auto& down_cans = subInfo.down_cans;
    c_idxs.resize(down_cans.size());
    p_idxs.resize(down_cans.size());
    c_offset.resize(down_cans.size()+1);
    p_offset.resize(down_cans.size()+1);
    down_group_num.clear();
    down_group_num.resize(down_cans.size(), 0);  // 0->no group
    for (ui i = 0; i < down_cans.size(); i ++) c_idxs[i] = i;
    c_cnt = 1;
    c_offset[0] = 0; c_offset[1] = down_cans.size();
    for (auto& unbr : down_indep) {
        vector<bool> down_grouped(c_idxs.size(), false);
        p_cnt = 0;
        for (ui i = 0; i < c_cnt; i++) {
            auto start = c_offset[i], end = c_offset[i+1];
            // 1. put all vertices that have no influence on unbr to one group
            p_cnt++;
            p_offset[p_cnt] = p_offset[p_cnt-1];
            for (auto idx1 = start; idx1 < end; idx1++) {
                auto& cans_idx = c_idxs[idx1];
                if (!influenced[unbr][cans_idx]) {
                    p_idxs[p_offset[p_cnt]++] = cans_idx;
                    down_grouped[cans_idx] = true;
                    down_group_num[cans_idx] = p_cnt;
                }
            }
            if (p_offset[p_cnt] == p_offset[p_cnt-1]) p_cnt--;
            for (auto idx1 = start; idx1 < end; idx1++) {
                auto& cans_idx1 = c_idxs[idx1];
                if (down_grouped[cans_idx1]) continue;
                p_cnt++;
                p_offset[p_cnt] = p_offset[p_cnt-1];
                p_idxs[p_offset[p_cnt]++] = cans_idx1;
                auto& cans1 = inf_cans[unbr][cans_idx1];
                down_grouped[cans_idx1] = true;
                down_group_num[cans_idx1] = p_cnt;
                for (auto idx2 = idx1 + 1; idx2 < end; idx2++) {
                    auto& cans_idx2 = c_idxs[idx2];
                    auto& cans2 = inf_cans[unbr][cans_idx2];
                    if (down_grouped[cans_idx2]) continue;
                    if (cans1.size() != cans2.size()) continue;
                    same = true;
                    for (ui com_idx = 0; com_idx < cans1.size(); com_idx++) {
                        if (cans1[com_idx] != cans2[com_idx]) {
                            same = false;
                            break;
                        }
                    }
                    if (same) {
                        p_idxs[p_offset[p_cnt]++] = cans_idx2;
                        down_grouped[cans_idx2] = true;
                        down_group_num[cans_idx2] = p_cnt;
                    }
                }
            }
        }
        c_idxs.swap(p_idxs);
        c_offset.swap(p_offset);
        c_cnt = p_cnt;
    }
    // after up_indep & down_indep info, compute edge info
    if (subInfo.connected) {
        comEdgeSubCon(index, depth);
    } else {
        comEdgeSubDis(index, depth);
    }
    return;
}

void
EvaluateQuery::comEdgeSubCon(FiPEIndex& index, ui depth) {
    auto& subInfo = index.subInfo_[depth];
    auto& subs = subInfo.subs;
    auto& up_cans = subInfo.up_cans;
    auto& up2down = subInfo.up2down;     // down_cans of each up_can
    auto& shared = subInfo.nbrs.shared_;
    auto& grouped = subInfo.grouped;
    auto& influenced = subInfo.influenced;
    auto& inf_cans = subInfo.inf_cans;
    auto& down_idxs = subInfo.down_idxs;
    auto& up_group_num = subInfo.subs.up_group_num;
    auto& down_group_num = subInfo.subs.down_group_num;

    // compute sub by use 2 offset alternately
    // c_* store the complete sub info, p_* store FiPE sub info
    auto& c_cnt = subs.c_cnt;
    auto& c_edge_up_offset = subs.edge_up_offset;
    auto& c_edge_up_idxs = subs.edge_up_idxs;
    auto& c_edge_down_offset = subs.edge_down_offset;
    auto& c_edge_down_idxs = subs.edge_down_idxs;
    auto& p_cnt = subs.p_cnt;
    auto& p_edge_up_offset = subs.p_edge_up_offset;
    auto& p_edge_up_idxs = subs.p_edge_up_idxs;
    auto& p_edge_down_offset = subs.p_edge_down_offset;
    auto& p_edge_down_idxs = subs.p_edge_down_idxs;
    ui max_size = up2down.size() > up_cans.size() ? up2down.size() : up_cans.size();

    // init edge_idxs
    c_cnt = 1;
    c_edge_up_offset.resize(max_size+2);
    p_edge_up_offset.resize(max_size+2);
    // c_edge_down_offset is initialized during edge building in comCurSpace
    c_edge_down_offset.resize(max_size+2);
    p_edge_down_offset.resize(max_size+2);
    p_edge_up_idxs.resize(max_size);
    c_edge_up_idxs.resize(max_size);
    c_edge_down_idxs.resize(up2down.size());
    p_edge_down_idxs.resize(up2down.size());
    c_edge_up_offset[0] = 0; c_edge_up_offset[1] = up_cans.size();
    for (ui i = 0; i < up_cans.size(); i++) c_edge_up_idxs[i] = i;
    for (ui i = 0; i < up2down.size(); i++) c_edge_down_idxs[i] = i;
    grouped.resize(up2down.size());
    bool same;

    for (auto& unbr : shared) {
        fill(grouped.begin(), grouped.end(), false);
        p_cnt = 0;
        for (ui group_idx = 0; group_idx < c_cnt; group_idx++) {
            // can not put all vertices that have no influence on unbr to one group
            // because down_group_num maybe different
            auto& up_start = c_edge_up_offset[group_idx];
            auto& up_end = c_edge_up_offset[group_idx+1];
            for (auto up_idx_idx1 = up_start; up_idx_idx1 < up_end; up_idx_idx1++) {
                auto& up_idx1 = c_edge_up_idxs[up_idx_idx1];
                auto& down_start1 = c_edge_down_offset[up_idx_idx1];
                auto& down_end1 = c_edge_down_offset[up_idx_idx1+1];
                for (auto down_idx_idx1 = down_start1; down_idx_idx1 < down_end1; down_idx_idx1++) {
                    auto& edge_idx1 = c_edge_down_idxs[down_idx_idx1];
                    if (grouped[edge_idx1]) continue;
                    auto& down_idx1 = down_idxs[edge_idx1];
                    auto& cans1 = inf_cans[unbr][edge_idx1];
                    p_cnt++;
                    p_edge_up_offset[p_cnt] = p_edge_up_offset[p_cnt-1];
                    p_edge_up_idxs[p_edge_up_offset[p_cnt]++] = up_idx1;
                    p_edge_down_offset[p_edge_up_offset[p_cnt]] = p_edge_down_offset[p_edge_up_offset[p_cnt] - 1];
                    p_edge_down_idxs[p_edge_down_offset[p_edge_up_offset[p_cnt]]++] = edge_idx1;
                    grouped[edge_idx1] = true;
                    for (ui up_idx_idx2 = up_idx_idx1; up_idx_idx2 < up_end; up_idx_idx2++) {
                        auto& up_idx2 = c_edge_up_idxs[up_idx_idx2];
                        auto& down_start2 = c_edge_down_offset[up_idx_idx2];
                        auto& down_end2 = c_edge_down_offset[up_idx_idx2+1];
                        if (up_idx2 != p_edge_up_idxs[p_edge_up_offset[p_cnt]-1]) {
                            p_edge_up_idxs[p_edge_up_offset[p_cnt]++] = up_idx2;
                            p_edge_down_offset[p_edge_up_offset[p_cnt]] = p_edge_down_offset[p_edge_up_offset[p_cnt] - 1];
                        }
                        for (auto down_idx_idx2 = down_start2; down_idx_idx2 < down_end2; down_idx_idx2++) {
                            auto& edge_idx2 = c_edge_down_idxs[down_idx_idx2];
                            auto& down_idx2 = down_idxs[edge_idx2];
                            if (grouped[edge_idx2]
                                || down_group_num[down_idx1] != down_group_num[down_idx2]
                                || up_group_num[up_idx1] != up_group_num[up_idx2]) continue;
                            same = true;
                            if (influenced[unbr][edge_idx2] || influenced[unbr][edge_idx1]) {
                                auto& cans2 = inf_cans[unbr][edge_idx2];
                                if (cans1.size() != cans2.size()) continue;
                                for (ui com_idx = 0; com_idx < cans1.size(); com_idx++) {
                                    if (cans1[com_idx] != cans2[com_idx]) {
                                        same = false;
                                        break;
                                    }
                                }
                            }
                            if (same) {
                                p_edge_down_idxs[p_edge_down_offset[p_edge_up_offset[p_cnt]]++] = edge_idx2;
                                grouped[edge_idx2] = true;
                            }
                        }
                        // if no edge is added for up_idx, edge_up_offset should not store this up_idx
                        if (p_edge_down_offset[p_edge_up_offset[p_cnt]] == p_edge_down_offset[p_edge_up_offset[p_cnt] - 1])
                            p_edge_up_offset[p_cnt]--;
                    }
                }
            }
        }
        c_edge_up_idxs.swap(p_edge_up_idxs);
        c_edge_up_offset.swap(p_edge_up_offset);
        c_edge_down_idxs.swap(p_edge_down_idxs);
        c_edge_down_offset.swap(p_edge_down_offset);
        c_cnt = p_cnt;
    }

    if (shared.empty()) {
        fill(grouped.begin(), grouped.end(), false);
        p_cnt = 0;
        for (ui group_idx = 0; group_idx < c_cnt; group_idx++) {
            auto& up_start = c_edge_up_offset[group_idx];
            auto& up_end = c_edge_up_offset[group_idx+1];
            for (auto up_idx_idx1 = up_start; up_idx_idx1 < up_end; up_idx_idx1++) {
                auto& up_idx1 = c_edge_up_idxs[up_idx_idx1];
                auto& down_start1 = c_edge_down_offset[up_idx_idx1];
                auto& down_end1 = c_edge_down_offset[up_idx_idx1+1];
                for (auto down_idx_idx1 = down_start1; down_idx_idx1 < down_end1; down_idx_idx1++) {
                    auto& edge_idx1 = c_edge_down_idxs[down_idx_idx1];
                    if (grouped[edge_idx1]) continue;
                    auto& down_idx1 = down_idxs[edge_idx1];
                    p_cnt++;
                    p_edge_up_offset[p_cnt] = p_edge_up_offset[p_cnt-1];
                    p_edge_up_idxs[p_edge_up_offset[p_cnt]++] = up_idx1;
                    p_edge_down_offset[p_edge_up_offset[p_cnt]] = p_edge_down_offset[p_edge_up_offset[p_cnt] - 1];
                    p_edge_down_idxs[p_edge_down_offset[p_edge_up_offset[p_cnt]]++] = edge_idx1;
                    grouped[edge_idx1] = true;
                    for (ui up_idx_idx2 = up_idx_idx1; up_idx_idx2 < up_end; up_idx_idx2++) {
                        auto& up_idx2 = c_edge_up_idxs[up_idx_idx2];
                        auto& down_start2 = c_edge_down_offset[up_idx_idx2];
                        auto& down_end2 = c_edge_down_offset[up_idx_idx2+1];
                        if (up_idx2 != p_edge_up_idxs[p_edge_up_offset[p_cnt]-1]) {
                            p_edge_up_idxs[p_edge_up_offset[p_cnt]++] = up_idx2;
                            p_edge_down_offset[p_edge_up_offset[p_cnt]] = p_edge_down_offset[p_edge_up_offset[p_cnt] - 1];
                        }
                        for (auto down_idx_idx2 = down_start2; down_idx_idx2 < down_end2; down_idx_idx2++) {
                            auto& edge_idx2 = c_edge_down_idxs[down_idx_idx2];
                            auto& down_idx2 = down_idxs[edge_idx2];
                            if (grouped[edge_idx2]
                                || down_group_num[down_idx1] != down_group_num[down_idx2]
                                || up_group_num[up_idx1] != up_group_num[up_idx2]) continue;
                                p_edge_down_idxs[p_edge_down_offset[p_edge_up_offset[p_cnt]]++] = edge_idx2;
                                grouped[edge_idx2] = true;
                        }
                        // if no edge is added for up_idx, edge_up_offset should not store this up_idx
                        if (p_edge_down_offset[p_edge_up_offset[p_cnt]] == p_edge_down_offset[p_edge_up_offset[p_cnt] - 1])
                            p_edge_up_offset[p_cnt]--;
                    }
                }
            }
        }
        c_edge_up_idxs.swap(p_edge_up_idxs);
        c_edge_up_offset.swap(p_edge_up_offset);
        c_edge_down_idxs.swap(p_edge_down_idxs);
        c_edge_down_offset.swap(p_edge_down_offset);
        c_cnt = p_cnt;
    }

    return;
}

void
EvaluateQuery::comEdgeSubDis(FiPEIndex& index, ui depth) {
    auto& subInfo = index.subInfo_[depth];
    auto& up_cans = subInfo.up_cans;
    auto& down_cans = subInfo.down_cans;
    auto& subs = subInfo.subs;
    auto& shared = subInfo.nbrs.shared_;
    auto& grouped = subInfo.grouped;
    auto& influenced = subInfo.influenced;
    auto& inf_cans = subInfo.inf_cans;
    auto& edge_valid = subInfo.edge_valid;
    auto& up_group_num = subInfo.subs.up_group_num;
    auto& down_group_num = subInfo.subs.down_group_num;

    // compute sub by use 2 offset alternately
    // c_* store the complete sub info, p_* store FiPE sub info
    auto& c_cnt = subs.c_cnt;
    auto& c_edge_up_offset = subs.edge_up_offset;
    auto& c_edge_up_idxs = subs.edge_up_idxs;
    auto& c_edge_down_offset = subs.edge_down_offset;
    auto& c_edge_down_idxs = subs.edge_down_idxs;
    auto& p_cnt = subs.p_cnt;
    auto& p_edge_down_offset = subs.p_edge_down_offset;
    auto& p_edge_down_idxs = subs.p_edge_down_idxs;

    // init edge_idxs
    auto up_num = up_cans.size();
    auto down_num = down_cans.size();
    auto edge_num = up_num*down_num;
    c_edge_down_offset.resize(edge_num+1);
    p_edge_down_offset.resize(edge_num+1);
    c_edge_down_idxs.resize(edge_num);
    p_edge_down_idxs.resize(edge_num);
    c_edge_down_offset[0] = 0, c_edge_down_offset[1] = edge_num;
    for (ui i = 0; i < edge_num; i++) c_edge_down_idxs[i] = i;
    grouped.resize(edge_num);
    bool same;

    // scanning edge, up_idx = edge_idx/down_cans.size(), down_idx = edge_idx%down_cans.size()
    //               up = up_cans[up_idx], down = down_cans[down_idx]
    for (auto& unbr : shared) {
        p_cnt = 0;
        fill(grouped.begin(), grouped.end(), false);
        for (ui group_id = 0; group_id < c_cnt; group_id++) {
            auto& edge_start = c_edge_down_offset[group_id];
            auto& edge_end = c_edge_down_offset[group_id+1];
            for (ui edge_idx_idx1 = edge_start; edge_idx_idx1 < edge_end; edge_idx_idx1++) {
                auto& edge_idx1 = c_edge_down_idxs[edge_idx_idx1];
                if (!edge_valid[edge_idx1] || grouped[edge_idx1]) continue;
                auto up_idx1 = edge_idx1/down_num;
                auto down_idx1 = edge_idx1%down_num;
                p_cnt+=1;
                p_edge_down_offset[p_cnt] = p_edge_down_offset[p_cnt - 1];
                p_edge_down_idxs[p_edge_down_offset[p_cnt]++] = edge_idx1;
                auto& cans1 = inf_cans[unbr][edge_idx1];
                grouped[edge_idx1] = true;
                for (ui edge_idx_idx2 = edge_idx_idx1 + 1; edge_idx_idx2 < edge_end; edge_idx_idx2++) {
                    auto& edge_idx2 = c_edge_down_idxs[edge_idx_idx2];
                    auto up_idx2 = edge_idx2/down_num;
                    auto down_idx2 = edge_idx2%down_num;
                    if (!edge_valid[edge_idx2] || grouped[edge_idx2]
                        || up_group_num[up_idx1] != up_group_num[up_idx2]
                        || down_group_num[down_idx1] != down_group_num[down_idx2]) continue;
                    same = true;
                    if (influenced[unbr][edge_idx2] || influenced[unbr][edge_idx1]) {
                        auto& cans2 = inf_cans[unbr][edge_idx2];
                        if (cans1.size() != cans2.size()) continue;
                        for (ui com_idx = 0; com_idx < cans1.size(); com_idx++) {
                            if (cans1[com_idx] != cans2[com_idx]) {
                                same = false;
                                break;
                            }
                        }
                    }
                    if (same) {
                        p_edge_down_idxs[p_edge_down_offset[p_cnt]++] = edge_idx2;
                        grouped[edge_idx2] = true;
                    }
                }
            }
        }
        c_edge_down_idxs.swap(p_edge_down_idxs);
        c_edge_down_offset.swap(p_edge_down_offset);
        c_cnt = p_cnt;
    }

    if (shared.empty()) {
        p_cnt = 0;
        fill(grouped.begin(), grouped.end(), false);
        for (ui group_id = 0; group_id < c_cnt; group_id++) {
            auto& edge_start = c_edge_down_offset[group_id];
            auto& edge_end = c_edge_down_offset[group_id+1];
            for (ui edge_idx_idx1 = edge_start; edge_idx_idx1 < edge_end; edge_idx_idx1++) {
                auto& edge_idx1 = c_edge_down_idxs[edge_idx_idx1];
                if (!edge_valid[edge_idx1] || grouped[edge_idx1]) continue;
                auto up_idx1 = edge_idx1/down_num;
                auto down_idx1 = edge_idx1%down_num;
                p_cnt+=1;
                p_edge_down_offset[p_cnt] = p_edge_down_offset[p_cnt - 1];
                p_edge_down_idxs[p_edge_down_offset[p_cnt]++] = edge_idx1;
                grouped[edge_idx1] = true;
                for (ui edge_idx_idx2 = edge_idx_idx1 + 1; edge_idx_idx2 < edge_end; edge_idx_idx2++) {
                    auto& edge_idx2 = c_edge_down_idxs[edge_idx_idx2];
                    auto up_idx2 = edge_idx2/down_num;
                    auto down_idx2 = edge_idx2%down_num;
                    if (!edge_valid[edge_idx2] || grouped[edge_idx2]
                        || up_group_num[up_idx1] != up_group_num[up_idx2]
                        || down_group_num[down_idx1] != down_group_num[down_idx2]) continue;
                        p_edge_down_idxs[p_edge_down_offset[p_cnt]++] = edge_idx2;
                        grouped[edge_idx2] = true;
                }
            }
        }
        c_edge_down_idxs.swap(p_edge_down_idxs);
        c_edge_down_offset.swap(p_edge_down_offset);
        c_cnt = p_cnt;
    }

    return;
}

/**
 * compute valid_cans of all unbrs for start_vertex
 */
bool
EvaluateQuery::comStartCans(FiPEIndex& index) {
    auto u = index.order_[0];
    auto& subInfo = index.subInfo_[0];
    auto& influenced = subInfo.influenced;
    auto& inf_cans = subInfo.inf_cans;
    auto& up_cans = subInfo.up_cans;
    auto& up_indep = subInfo.nbrs.up_indep_;
    auto& cans = index.valid_cans_[u].back();
    auto can_cnt = index.valid_cans_[u].back().size();
    for (auto& unbr : up_indep) {
        influenced[unbr].clear();
        inf_cans[unbr].clear();
        influenced[unbr].reserve(can_cnt);
        inf_cans[unbr].reserve(can_cnt);
    }
    up_cans.clear();
    up_cans.reserve(can_cnt);

    // compute the valid_cans of all nbrs for each cans
    ui unbrs_cnt;
    auto unbrs = index.q_graph_->getVertexNeighbors(u, unbrs_cnt);
    for (ui v_idx = 0; v_idx < can_cnt; v_idx++) {
        auto& v = cans[v_idx];
        bool valid = true;
        ui edge_idx = v_idx;
        if (index.subCans_[0].splitted) edge_idx += index.subCans_[0].cur_start;
        for (ui i = 0; i < unbrs_cnt; i++) {
            auto& unbr = unbrs[i];
            // old valid_cans of unbr
            auto& edges = *(index.index_[u][unbr].back());
            auto vnbrs_cnt = edges.offset_[edge_idx + 1] - edges.offset_[edge_idx];
            if (vnbrs_cnt == 0) {
                valid = false;
                break;
            }
        }
        if (valid) {
            up_cans.emplace_back(cans[v_idx]);
            for (auto& unbr : up_indep) {
                auto& edges = *(index.index_[u][unbr].back());
                auto vnbrs_cnt = edges.offset_[edge_idx + 1] - edges.offset_[edge_idx];
                auto vnbrs = edges.edge_ + edges.offset_[edge_idx];
                influenced[unbr].emplace_back(true);  // can be removed for start_vertex
                inf_cans[unbr].emplace_back(vnbrs, vnbrs+vnbrs_cnt);
            }
        }
    }
    return up_cans.size();  // 0->false, o.w.->true
}

/**
 * compute valid_cans for down_indep_ & shared_
 */
bool
EvaluateQuery::comCurSpace(FiPEIndex& index, ui depth) {
    auto& up = index.order_[depth];
    auto& down = index.order_[depth+1];
    auto& subInfo = index.subInfo_[depth];
    auto& up_cans = subInfo.up_cans;
    auto& subs = subInfo.subs;
    auto& edge_down_offset = subs.edge_down_offset;  // valid_cans range of each cans to nxt
    auto& up2down = subInfo.up2down;       // valid_cans of nxt
    auto& down_record = subInfo.subs.down_record;
    bool success = true;

    // 1.compute the valid_cans of shared for up
    auto& up_shared_valid_cans = subInfo.up_shared_valid_cans;      // cans of [unbr][up_idx]
    auto& shared = subInfo.nbrs.shared_;
    ui up_valid_idx = 0;
    for (auto& nbr : shared) {
        up_shared_valid_cans[nbr].resize(up_cans.size());
    }
    for (ui v_idx = 0; v_idx < up_cans.size(); v_idx++, up_valid_idx++) {
        up_cans[up_valid_idx] = up_cans[v_idx];
        for (auto& unbr : shared) {
            auto vnbrs = move(index.getNeighbors(up, unbr, up_cans[v_idx]));
            if (vnbrs.size() == 0) {
                up_valid_idx--;
                break;
            }
            up_shared_valid_cans[unbr][up_valid_idx].swap(vnbrs);
        }
    }
    if (up_valid_idx == 0) return false;
    up_cans.resize(up_valid_idx);

    // build edges
    if (subInfo.connected) {
        edge_down_offset.clear();
        edge_down_offset.emplace_back(0);
        up2down.clear();
        for (auto& up_can : up_cans) {
            auto vnbrs = move(index.getNeighbors(up, down, up_can));
            edge_down_offset.emplace_back(edge_down_offset.back()+vnbrs.size());
            up2down.insert(up2down.end(), move(vnbrs.begin()), move(vnbrs.end()));
        }
    }

    // compute the down_cans
    auto& down_cans = subInfo.down_cans;
    if (subInfo.connected) {
        down_cans.clear();
        auto& down_idxs = subInfo.down_idxs;
        down_idxs.resize(up2down.size());
        for (ui i = 0; i < up2down.size(); i++) {
            if (down_record[up2down[i]] == (ui)-1) {
                down_record[up2down[i]] = down_cans.size();
                down_cans.emplace_back(up2down[i]);
            }
            down_idxs[i] = down_record[up2down[i]];
        }
    } else {
        down_cans = index.valid_cans_[down].back();
    }

    // compute valid_cans for down_indep & down_shared_valid_cans & delayed_nbrs
    auto& influenced = subInfo.influenced;
    auto& inf_cans = subInfo.inf_cans;
    auto& down_indep = subInfo.nbrs.down_indep_;
    auto& delayed_nbrs = subInfo.nbrs.delayed_;
    ui down_valid_idx = 0;
    auto& down_valid = subInfo.down_valid;
    auto& down_shared_valid_cans = subInfo.down_shared_valid_cans;
    down_valid.clear();
    down_valid.resize(down_cans.size(), true);
    for (auto& nbr : down_indep) {
        influenced[nbr].clear();
        inf_cans[nbr].clear();
        influenced[nbr].resize(down_cans.size());
        inf_cans[nbr].resize(down_cans.size());
    }
    for (auto& nbr : shared) {
        down_shared_valid_cans[nbr].resize(down_cans.size());
    }
    for (auto& nbr : delayed_nbrs) {
        influenced[nbr].clear();
        inf_cans[nbr].clear();
        influenced[nbr].resize(down_cans.size());
        inf_cans[nbr].resize(down_cans.size());
    }
    for (ui v_idx = 0; v_idx < down_cans.size(); v_idx++, down_valid_idx++) {
        auto& v = down_cans[v_idx];
        down_record[v] = down_valid_idx;
        for (auto& unbr : down_indep) {
            auto vnbrs = move(index.getNeighbors(down, unbr, v));
            if (vnbrs.size() == 0) {
                down_valid[v_idx] = false;
                down_valid_idx--;
                break;
            }
            if (index.valid_cans_[unbr].back().size() == vnbrs.size()) {
                influenced[unbr][down_valid_idx] = false;
                continue;  // means no changes
            }

            // write to down_valid_idx, down_cans will be shrinked at the end of this function
            influenced[unbr][down_valid_idx] = true;
            inf_cans[unbr][down_valid_idx].swap(vnbrs);
        }
        if (!down_valid[v_idx]) continue;
        for (auto& unbr : shared) {
            auto vnbrs = move(index.getNeighbors(down, unbr, v));
            if (vnbrs.size() == 0) {
                down_valid_idx--;
                down_valid[v_idx] = false;
                break;
            }
            down_shared_valid_cans[unbr][down_valid_idx].swap(vnbrs);
        }
        if (!down_valid[v_idx]) continue;
        for (auto& unbr : delayed_nbrs) {
            auto vnbrs = move(index.getNeighbors(down, unbr, v));
            if (vnbrs.size() == 0) {
                down_valid[v_idx] = false;
                down_valid_idx--;
                break;
            }
            if (index.valid_cans_[unbr].back().size() == vnbrs.size()) {
                influenced[unbr][down_valid_idx] = false;
                continue;  // means no changes
            }

            influenced[unbr][down_valid_idx] = true;
            inf_cans[unbr][down_valid_idx].swap(vnbrs);
        }
    }
    if (down_valid_idx != 0) {
        // compute the valid_cans of shared for edges, based on down(up)_shared_valid_cans
        if (subInfo.connected) {
            success = comSharedCon(index, depth);
        } else {
            success = comSharedDis(index, depth);
        }
    } else {
        success = false;
    }

    // recover down_record
    ui valid_idx = 0;
    for (ui v_idx = 0; v_idx < down_cans.size(); v_idx++) {
        auto& down_can = down_cans[v_idx];
        down_record[down_can] = (ui)-1;
        if (down_valid[v_idx]) down_cans[valid_idx++] = down_can;
    }
    down_cans.resize(valid_idx);
    return success;
}

inline bool
EvaluateQuery::comSharedCon(FiPEIndex& index, ui depth) {
    auto& subInfo = index.subInfo_[depth];
    auto& up_cans = subInfo.up_cans;
    auto& subs = subInfo.subs;
    auto& edge_down_offset = subs.edge_down_offset;  // valid_cans range of each cans to nxt
    auto& up2down = subInfo.up2down;                // valid_cans of nxt
    // set group[idx] to true to jump the grouping of invalid edges
    auto& shared = subInfo.nbrs.shared_;
    auto& edge_valid = subInfo.edge_valid;
    auto& influenced = subInfo.influenced;
    auto& inf_cans = subInfo.inf_cans;
    edge_valid.clear();
    edge_valid.resize(up2down.size(), true);
    auto& down_valid = subInfo.down_valid;
    auto& down_record = subInfo.subs.down_record;  // record v_idx to valid_down_cans
    auto& up_shared_valid_cans = subInfo.up_shared_valid_cans;      // cans of [unbr][up_idx]
    auto& down_shared_valid_cans = subInfo.down_shared_valid_cans;  // cans of [unbr][down_idx]
    auto& down_idxs = subInfo.down_idxs;
    auto& up_idxs = subInfo.up_idxs;
    up_idxs.resize(down_idxs.size());

    // reserve spaces for influenced
    for (auto& nbr : shared) {
        influenced[nbr].clear();
        inf_cans[nbr].clear();
        influenced[nbr].resize(up2down.size());
        inf_cans[nbr].resize(up2down.size());
    }

    // compute shared_valid_cans,
    // by the way, remove invalid down_cans & edges, update edge_down_offset, up2down, down_idxs
    ui edge_valid_idx = 0;
    auto valid_edge_down_offset = edge_down_offset;
    for (ui up_idx = 0; up_idx < up_cans.size(); up_idx++) {
        auto& edge_start = edge_down_offset[up_idx];
        auto& edge_end = edge_down_offset[up_idx+1];
        valid_edge_down_offset[up_idx+1] = valid_edge_down_offset[up_idx];
        for (ui edge_idx = edge_start; edge_idx < edge_end; edge_idx++) {
            auto& down = up2down[edge_idx];
            auto& valid_down_can_idx = down_record[down];
            auto& down_idx = down_idxs[edge_idx];
            if (!down_valid[down_idx]) continue;
            for (auto& nbr : shared) {
                auto vnbrs = move(SetOp::intersectTwo(up_shared_valid_cans[nbr][up_idx],
                                                      down_shared_valid_cans[nbr][valid_down_can_idx]));
                if (vnbrs.size() == 0) {
                    edge_valid[edge_idx] = false;
                    break;
                }
                if (vnbrs.size() == index.valid_cans_[nbr].back().size()) {
                    influenced[nbr][edge_valid_idx] = false;
                    continue;
                }
                influenced[nbr][edge_valid_idx] = true;
                inf_cans[nbr][edge_valid_idx].swap(vnbrs);
            }
            if (edge_valid[edge_idx]) {
                edge_valid_idx++;
                up2down[valid_edge_down_offset[up_idx+1]] = down;
                down_idxs[valid_edge_down_offset[up_idx+1]] = valid_down_can_idx;
                up_idxs[valid_edge_down_offset[up_idx+1]] = up_idx;
                valid_edge_down_offset[up_idx+1]++;
            }
        }
    }
    edge_down_offset.swap(valid_edge_down_offset);
    up2down.resize(edge_valid_idx);
    return up2down.size();
}

inline bool
EvaluateQuery::comSharedDis(FiPEIndex& index, ui depth) {
    auto& subInfo = index.subInfo_[depth];
    auto& up_cans = subInfo.up_cans;
    auto& down_cans = subInfo.down_cans;
    // set group[idx] to true to jump the grouping of invalid edges
    auto& shared = subInfo.nbrs.shared_;
    auto& edge_valid = subInfo.edge_valid;
    auto& influenced = subInfo.influenced;
    auto& inf_cans = subInfo.inf_cans;
    auto& down_valid = subInfo.down_valid;
    auto& up_valid_cans = subInfo.up_shared_valid_cans;      // cans of [unbr][up_idx]
    auto& down_valid_cans = subInfo.down_shared_valid_cans;  // cans of [unbr][down_idx]
    auto valid_down_cans = down_cans;
    ui valid_idx = 0;
    for (ui v_idx = 0; v_idx < down_cans.size(); v_idx++) {
        auto& down_can = down_cans[v_idx];
        if (down_valid[v_idx]) valid_down_cans[valid_idx++] = down_can;
    }
    valid_down_cans.resize(valid_idx);
    auto edge_num = up_cans.size()*valid_down_cans.size();
    edge_valid.clear();
    edge_valid.resize(edge_num, true);
    for (auto& nbr : shared) {
        influenced[nbr].clear();
        inf_cans[nbr].clear();
        influenced[nbr].resize(edge_num);
        inf_cans[nbr].resize(edge_num);
    }

    for (ui up_idx = 0; up_idx < up_cans.size(); up_idx++) {
        for (ui down_idx = 0; down_idx < valid_down_cans.size(); down_idx++) {
            auto edge_idx = up_idx*valid_down_cans.size() + down_idx;
            for (auto& nbr : shared) {
                auto vnbrs = move(SetOp::intersectTwo(up_valid_cans[nbr][up_idx], down_valid_cans[nbr][down_idx]));
                if (vnbrs.size() == 0) {
                    edge_valid[edge_idx] = false;
                    break;
                }
                if (vnbrs.size() == index.valid_cans_[nbr].back().size()) {
                    influenced[nbr][edge_idx] = false;
                    continue;
                }
                influenced[nbr][edge_idx] = true;
                inf_cans[nbr][edge_idx].swap(vnbrs);
            }
        }
    }
    return true;
}
