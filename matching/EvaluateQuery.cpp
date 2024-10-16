#include "EvaluateQuery.h"
#include "utility/bsx/IndepSet.h"
#include "utility/bsx/nodeSim.h"
#include "utility/bsx/SetOp.h"
#include <stack>
#include <vector>
#include <cstring>
#include <sys/stat.h>

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

void
EvaluateQuery::BS1Engine(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix,
                             ui **candidates, ui *candidates_count, ui *order, ui *pivot,
                             size_t output_limit_num, size_t &call_count, mpz_t embedding_cnt) {
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
    mpz_init_set_ui(embedding_cnt, 0);
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

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
    out_files.resize(max_depth);
    for (ui i = 1; i < max_depth; i++) {
        out_files[i-1].open(g_name + "/" + std::to_string(i-1) + ".txt");
    }
    memset(embedding, (ui)-1, sizeof(ui)*max_depth);
    std::vector<std::stack<std::vector<ui>>> valid_cans;
    for (ui i = 0; i < max_depth; i++) {
        std::stack<std::vector<ui>> sv;
        std::vector<ui> vec;
        vec.insert(vec.end(), candidates[i], candidates[i]+candidates_count[i]);
        sv.push(std::move(vec));
        valid_cans.emplace_back(std::move(sv));
    }
#endif

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

            if (cur_depth == max_depth - 1) {
                mpz_add_ui(embedding_cnt, embedding_cnt, 1);
                visited_vertices[v] = false;
                if (output_limit_num != (size_t)-1 && mpz_cmp_ui(embedding_cnt, output_limit_num) > 0) {
                    goto EXIT;
                }
            } else {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidateIndex(data_graph, cur_depth, embedding, idx_embedding, idx_count,
                                            valid_candidate_idx, edge_matrix, visited_vertices, bn,
                                            bn_count, order, pivot, candidates, query_graph);
#ifdef ANALYZE_DUPLICATE
                // compute valid_cans of all u's connected to u
                ui unbrs_cnt = 0;
                auto unbrs = query_graph->getVertexNeighbors(u, unbrs_cnt);
                for (ui i = 0; i < unbrs_cnt; i++) {
                    auto unbr = unbrs[i];
                    if (embedding[unbr] != (ui)-1) {  // if not matched, compute valid_cans
                        // first, get all neighbors of v on dataGraph
                        ui vnbrs_cnt = 0;
                        auto vnbrs = data_graph->getVertexNeighbors(v, vnbrs_cnt);
                        // and then intersected with valid_cans[unbr] & push
                        auto res = SetOp::intersectTwo(valid_cans[unbr].top(), vnbrs, vnbrs_cnt);
                        valid_cans[unbr].push(std::move(res));
                    }
                }
                for (ui i = cur_depth; i < max_depth; i++) {
                    auto cur_u = order[i];
                    for (auto& can : valid_cans[i].top()) {
                        out_files[cur_depth-1] << can << " ";
                    }
                    out_files[cur_depth-1] << std::endl;
                }
                out_files[cur_depth-1] << "------" << std::endl;
#endif
            }
        }

        // backtrack
        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
#ifdef ANALYZE_DUPLICATE
        // if (idx_count[cur_depth+1] == 0) continue;
        auto last_u = order[cur_depth + 1];
        embedding[last_u] = (ui)-1;
        auto u = order[cur_depth];
        auto v = embedding[u];
        // restore of neighbors of u valid_cans
        ui unbrs_cnt = 0;
        auto unbrs = query_graph->getVertexNeighbors(u, unbrs_cnt);
        for (ui i = 0; i < unbrs_cnt; i++) {
            auto unbr = unbrs[i];
            if (embedding[unbr] != (ui)-1) {
                valid_cans[unbr].pop();
            }
        }
#endif
    }


    // Release the buffer.
    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);
#ifdef ANALYZE_DUPLICATE
    for (ui i = 1; i < max_depth; i++) { 
        out_files[i-1].close();
    }
#endif

    return;
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
                    if (!data_graph->checkEdgeExistence(temp_v, u_bn_v)) {
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

/**
 * use bsx method
*/
void
EvaluateQuery::BSXEngine(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix,
                          ui **candidates, ui *candidates_count, ui *order,
                          size_t output_limit_num, size_t &call_count, mpz_t embedding_cnt) {
    ui q_num = query_graph->getVerticesCount();
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
    bsxComEqBatch(index, start_vertex);
    batch_info[start_vertex].print();
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
            ui cur_batch_cnt;
            VertexID* cur_batch = batch_info[u].cur_batch(cur_batch_cnt);

            // nxt batch
            batch_info[depth2u[cur_depth]].idx_.top()++;

            auto& cur_cans_cnt = index.valid_cnt_[u].top();
            auto& cur_cans = index.valid_cans_[u].top();
            std::copy(cur_batch, cur_batch+cur_batch_cnt, cur_cans);
            cur_cans_cnt = cur_batch_cnt;

            VertexID failed_u = bsxRefine(index, u);
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
                bsxGenResult(num_indep, indep_nodes, index);
                gmp_printf("new result: %Zd\n", level_embeddings);
                mpz_add(embedding_cnt, embedding_cnt, level_embeddings);
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
                bsxComEqBatch(index, cur_u);
                batch_info[cur_u].print();
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
        VertexID cur_u = depth2u[cur_depth];
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

    return;
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
    std::copy(index.visited_u, index.visited_u+q_num, nbr_updated);

    // process first u in influenced_u seperately, valid_cans of first_influenced_u comes from
    //   its batch_nodes, and couldn't be deleted
    auto first_u = index.influenced_u_.top()[0];
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
        batch_idx[i] = batch_cnt++;
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
                batches[offset[num]+(cnt[num]++)] = nodes[ev_idx];
            }
        }
        // process max_cnt
        if (maxCnt < cnt[num]) maxCnt = cnt[num];
        // no need for sort, because the idxs have been sorted and this fun. will not break it
        // std::sort(batch+offset[num], batch+(offset[num]+cnt[num]));
        num++;
    }
    delete[] batch_idx;
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
    index.updateIndex(influenced, u);
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
