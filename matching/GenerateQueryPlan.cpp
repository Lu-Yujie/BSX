#include "GenerateQueryPlan.h"
#include "computesetintersection.h"
#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>
#include "graphoperations.h"
#include <memory.h>

void GenerateQueryPlan::generateGQLQueryPlan(const Graph *data_graph, const Graph *query_graph, ui *candidates_count,
                                             ui *&order, ui *&pivot) {
     /**
      * Select the vertex v such that (1) v is adjacent to the selected vertices; and (2) v has the minimum number of candidates.
      */
     std::vector<bool> visited_vertices(query_graph->getVerticesCount(), false);
     std::vector<bool> adjacent_vertices(query_graph->getVerticesCount(), false);
     order = new ui[query_graph->getVerticesCount()];
     pivot = new ui[query_graph->getVerticesCount()];

     VertexID start_vertex = selectGQLStartVertex(query_graph, candidates_count);
     order[0] = start_vertex;
     updateValidVertices(query_graph, start_vertex, visited_vertices, adjacent_vertices);

     for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
          VertexID next_vertex = 0;
          ui min_value = data_graph->getVerticesCount() + 1;
          for (ui j = 0; j < query_graph->getVerticesCount(); ++j) {
               VertexID cur_vertex = j;

               if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex]) {
                    if (candidates_count[cur_vertex] < min_value) {
                         min_value = candidates_count[cur_vertex];
                         next_vertex = cur_vertex;
                    }
                    else if (candidates_count[cur_vertex] == min_value && query_graph->getVertexDegree(cur_vertex) > query_graph->getVertexDegree(next_vertex)) {
                         next_vertex = cur_vertex;
                    }
               }
          }
          updateValidVertices(query_graph, next_vertex, visited_vertices, adjacent_vertices);
          order[i] = next_vertex;
     }

     // Pick a pivot randomly.
     for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
         VertexID u = order[i];
         for (ui j = 0; j < i; ++j) {
             VertexID cur_vertex = order[j];
             if (query_graph->checkEdgeExistence(u, cur_vertex)) {
                 pivot[i] = cur_vertex;
                 break;
             }
         }
     }
}

VertexID GenerateQueryPlan::selectGQLStartVertex(const Graph *query_graph, ui *candidates_count) {
    /**
     * Select the vertex with the minimum number of candidates as the start vertex.
     * Tie Handling:
     *  1. degree
     *  2. label id
     */

     ui start_vertex = 0;

     for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
          VertexID cur_vertex = i;

          if (candidates_count[cur_vertex] < candidates_count[start_vertex]) {
               start_vertex = cur_vertex;
          }
          else if (candidates_count[cur_vertex] == candidates_count[start_vertex]
                   && query_graph->getVertexDegree(cur_vertex) > query_graph->getVertexDegree(start_vertex)) {
               start_vertex = cur_vertex;
          }
     }

     return start_vertex;
}

void GenerateQueryPlan::updateValidVertices(const Graph *query_graph, VertexID query_vertex, std::vector<bool> &visited,
                                            std::vector<bool> &adjacent) {
     visited[query_vertex] = true;
     ui nbr_cnt;
     const ui* nbrs = query_graph->getVertexNeighbors(query_vertex, nbr_cnt);

     for (ui i = 0; i < nbr_cnt; ++i) {
          ui nbr = nbrs[i];
          adjacent[nbr] = true;
     }
}

void GenerateQueryPlan::checkQueryPlanCorrectness(const Graph *query_graph, ui *order, ui *pivot) {
    ui q_num = query_graph->getVerticesCount();
    std::vector<bool> visited_vertices(q_num, false);
    // Check whether each query vertex is in the order.
    for (ui i = 0; i < q_num; ++i) {
        VertexID vertex = order[i];
        assert(vertex < q_num && vertex >= 0);

        visited_vertices[vertex] = true;
    }

    for (ui i = 0; i < q_num; ++i) {
        VertexID vertex = i;
        assert(visited_vertices[vertex]);
    }

    // Check whether the order is connected.

    std::fill(visited_vertices.begin(), visited_vertices.end(), false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < q_num; ++i) {
        VertexID vertex = order[i];
        VertexID pivot_vertex = pivot[i];
        assert(query_graph->checkEdgeExistence(vertex, pivot_vertex));
        assert(visited_vertices[pivot_vertex]);
        visited_vertices[vertex] = true;
    }
}

void GenerateQueryPlan::printQueryPlan(const Graph *query_graph, ui *order) {
    ui q_num = query_graph->getVerticesCount();
    printf("Query Plan: ");
    for (ui i = 0; i < q_num; ++i) {
        printf("%u, ", order[i]);
    }
    printf("\n");

    printf("%u: N/A\n", order[0]);
    for (ui i = 1; i < q_num; ++i) {
        VertexID end_vertex = order[i];
        printf("%u: ", end_vertex);
        for (ui j = 0; j < i; ++j) {
            VertexID begin_vertex = order[j];
            if (query_graph->checkEdgeExistence(begin_vertex, end_vertex)) {
                printf("R(%u, %u), ", begin_vertex, end_vertex);
            }
        }
        printf("\n");
    }
}

void GenerateQueryPlan::printSimplifiedQueryPlan(const Graph *query_graph, ui *order) {
    ui q_num = query_graph->getVerticesCount();
    printf("Query Plan: ");
    for (ui i = 0; i < q_num; ++i) {
        printf("%u ", order[i]);
    }
    printf("\n");
}

void GenerateQueryPlan::generateDSPisoQueryPlan(const Graph *query_graph, Edges ***edge_matrix, ui *&order, ui *&pivot,
                                                TreeNode *tree, ui *bfs_order, ui *candidates_count, ui **&weight_array) {
    ui q_num = query_graph->getVerticesCount();
    order = new ui[q_num];
    pivot = new ui[q_num];

    for (ui i = 0; i < q_num; ++i) {
        order[i] = bfs_order[i];
    }

    for (ui i = 1; i < q_num; ++i) {
        pivot[i] = tree[order[i]].parent_;
    }

    // Compute weight array.
    weight_array = new ui*[q_num];
    for (ui i = 0; i < q_num; ++i) {
        weight_array[i] = new ui[candidates_count[i]];
        std::fill(weight_array[i], weight_array[i] + candidates_count[i], std::numeric_limits<ui>::max());
    }

    for (int i = q_num - 1; i >= 0; --i) {
        VertexID vertex = order[i];
        TreeNode& node = tree[vertex];
        bool set_to_one = true;

        for (ui j = 0; j < node.fn_count_; ++j) {
            VertexID child = node.fn_[j];
            TreeNode& child_node = tree[child];

            if (child_node.bn_count_ == 1) { // the tree-like contraint
                set_to_one = false;
                Edges& cur_edge = *edge_matrix[vertex][child];
                for (ui k = 0; k < candidates_count[vertex]; ++k) {
                    ui cur_candidates_count = cur_edge.offset_[k + 1] - cur_edge.offset_[k];
                    ui* cur_candidates = cur_edge.edge_ + cur_edge.offset_[k];

                    ui weight = 0;

                    for (ui l = 0; l < cur_candidates_count; ++l) {
                        ui candidates = cur_candidates[l];
                        weight += weight_array[child][candidates];
                    }

                    if (weight < weight_array[vertex][k])
                        weight_array[vertex][k] = weight;
                }
            }
        }

        if (set_to_one) {
            std::fill(weight_array[vertex], weight_array[vertex] + candidates_count[vertex], 1);
        }
    }
}

void GenerateQueryPlan::checkQueryPlanCorrectness(const Graph *query_graph, ui *order) {
    ui q_num = query_graph->getVerticesCount();
    std::vector<bool> visited_vertices(q_num, false);
    // Check whether each query vertex is in the order.
    for (ui i = 0; i < q_num; ++i) {
        VertexID vertex = order[i];
        assert(vertex < q_num && vertex >= 0);

        visited_vertices[vertex] = true;
    }

    for (ui i = 0; i < q_num; ++i) {
        VertexID vertex = i;
        assert(visited_vertices[vertex]);
    }

    // Check whether the order is connected.

    std::fill(visited_vertices.begin(), visited_vertices.end(), false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < q_num; ++i) {
        VertexID u = order[i];

        bool valid = false;
        for (ui j = 0; j < i; ++j) {
            VertexID v = order[j];
            if (query_graph->checkEdgeExistence(u, v)) {
                valid = true;
                break;
            }
        }

        assert(valid);
        visited_vertices[u] = true;
    }
}

void
GenerateQueryPlan::generateRMQueryPlan(const Graph *query_graph, ui *&order, Edges ***edge_matrix, ui *&pivot) {
    std::vector<nd_tree_node> density_tree;
    std::vector<nd_tree_node> k12_tree;
    std::vector<nd_tree_node> k23_tree;
    std::vector<nd_tree_node> k34_tree;

    nd_interface::nd(query_graph, 1, 2, k12_tree);
    nd_interface::nd(query_graph, 2, 3, k23_tree);
    nd_interface::nd(query_graph, 3, 4, k34_tree);

    construct_density_tree(density_tree, k12_tree, k23_tree, k34_tree);

    std::vector<std::vector<uint32_t>> node_orders;
    std::vector<std::vector<uint32_t>> vertex_orders;
    traversal_density_tree(query_graph, edge_matrix, density_tree, vertex_orders, node_orders);
    
    select_vertex_order(query_graph, vertex_orders, node_orders);
    if (order == NULL) {
        order = new ui[vertex_orders.back().size()];
    }
    if (pivot == NULL) {
        pivot = new ui[vertex_orders.back().size()];
    }
    memcpy(order, vertex_orders.back().data(), sizeof(ui)*vertex_orders.back().size());
    // compute pivot
    for (ui i = 1; i < query_graph->getVerticesCount(); i++) {
        for (ui j = 0; j < i; ++j) {
            if (query_graph->checkEdgeExistence(order[i], order[j])) {
                pivot[i] = order[j];
                break;
            }
        }
    }
}

void GenerateQueryPlan::construct_density_tree(std::vector<nd_tree_node> &density_tree,
                                               std::vector<nd_tree_node> &k12_tree,
                                               std::vector<nd_tree_node> &k23_tree,
                                               std::vector<nd_tree_node> &k34_tree) {
    /** Eliminate the duplicate nodes.
     *     1. For each nucleus in k12_tree, if it is a subset of the nucleus in k23, then remove the subtree rooted at it from k12_tree.
     *     2. For each nucleus in k23_tree, if it is a subset of the nucleus in k34, then remove the subtree rooted at it from k23_tree.
     */
    eliminate_node(density_tree, k12_tree, k23_tree);
    eliminate_node(density_tree, k23_tree, k34_tree);

    /**
     * Add the node in k34_tree to the density tree.
     */
    uint32_t offset = density_tree.size();
    for (auto& node : k34_tree) {
        density_tree.push_back(node);
        auto& new_node = density_tree.back();

        if (new_node.parent_ != -1) {
            new_node.parent_ += offset;
        }

        new_node.id_ += offset;

        for (uint32_t i = 0; i < new_node.children_.size(); ++i) {
            new_node.children_[i] += offset;
        }
    }

    /**
     * Construct the tree.
     *     1. For the root node of each k23_tree, find the smallest (by cardinality) 12 nucleus contain it and set the nucleus as its parent.
     *     2. For the root node of each k34_tree, find the smallest (by cardinality) 12 or 23 nucleus contain it and set the nucleus as its parent.
     */
     merge_tree(density_tree);
}

void GenerateQueryPlan::eliminate_node(std::vector<nd_tree_node> &density_tree,
                                       std::vector<nd_tree_node> &src_tree,
                                       std::vector<nd_tree_node> &target_tree) {
    std::vector<int> roots;
    for (auto& node : target_tree) {
        if (node.parent_ == -1) {
            roots.push_back(node.id_);
        }
    }

    std::queue<int> q;
    unordered_map<int, int> mappings;
    mappings[-1] = -1;

    for (auto& node : src_tree) {
        if (node.parent_ == -1) {
            q.push(node.id_);
            while (!q.empty()) {
                int cur_node_id = q.front();
                q.pop();

                auto &cur_node = src_tree[cur_node_id];

                bool remove = false;
                for (auto &target_node_id : roots) {
                    auto &target_node = target_tree[target_node_id];
                    if (cur_node.vertices_.size() <= target_node.vertices_.size()) {
                        uint32_t cn_cnt = 0;
                        ComputeSetIntersection::ComputeCandidates((uint32_t *) cur_node.vertices_.data(),
                                                                  (uint32_t) cur_node.vertices_.size(),
                                                                  (uint32_t *) target_node.vertices_.data(),
                                                                  (uint32_t) target_node.vertices_.size(), cn_cnt);

                        if (cn_cnt == cur_node.vertices_.size()) {
                            remove = true;
                        }
                    }
                }

                if (!remove) {
                    auto new_id = static_cast<int>(density_tree.size());
                    mappings[cur_node.id_] = new_id;
                    density_tree.push_back(cur_node);
                    // Update its id, parent id and the children of its parent.
                    density_tree.back().id_ = new_id;
                    density_tree.back().parent_ = mappings[cur_node.parent_];
                    density_tree.back().children_.clear();

                    if (density_tree.back().parent_ != -1)
                        density_tree[density_tree.back().parent_].children_.push_back(new_id);

                    for (auto child_id : cur_node.children_) {
                        q.push(child_id);
                    }
                }
            }
        }
    }
}

void GenerateQueryPlan::merge_tree(std::vector<nd_tree_node> &density_tree) {
    for (int i = density_tree.size() - 1; i >= 0; --i) {
        auto& cur_node = density_tree[i];
        uint32_t smallest_cardinality = std::numeric_limits<uint32_t>::max();
        if (cur_node.parent_ == -1) {
            for (int j = i - 1; j >= 0; --j) {
                auto& target_node = density_tree[j];

                if (target_node.r_ < cur_node.r_ && target_node.vertices_.size() < smallest_cardinality) {
                    uint32_t cn_cnt = 0;
                    ComputeSetIntersection::ComputeCandidates((uint32_t *) cur_node.vertices_.data(),
                                                              (uint32_t) cur_node.vertices_.size(),
                                                              (uint32_t *) target_node.vertices_.data(),
                                                              (uint32_t) target_node.vertices_.size(), cn_cnt);

                    if (cn_cnt == cur_node.vertices_.size()) {
                        cur_node.parent_ = target_node.id_;
                        target_node.children_.push_back(cur_node.id_);
                        smallest_cardinality = target_node.vertices_.size();
                    }
                }

            }
        }
    }
}

void GenerateQueryPlan::traversal_density_tree(const Graph *query_graph, Edges***edge_matrix,
                                               std::vector<nd_tree_node> &density_tree,
                                               std::vector<std::vector<uint32_t>> &vertex_orders,
                                               std::vector<std::vector<uint32_t>> &node_orders) {
    std::vector<uint32_t> vertex_order;
    std::vector<uint32_t> node_order;
    std::vector<bool> visited_vertex(query_graph->getVerticesCount(), false);
    std::vector<bool> visited_node(query_graph->getVerticesCount(), false);
    std::unordered_set<uint32_t> extendable_vertex;

    for (auto& node : density_tree) {
        if (node.children_.empty()) {
            traversal_node(query_graph, edge_matrix, density_tree, vertex_order, node_order, visited_vertex, visited_node,
                           node, extendable_vertex);

            vertex_orders.emplace_back(vertex_order);
            node_orders.emplace_back(node_order);

            vertex_order.clear();
            node_order.clear();
            extendable_vertex.clear();

            std::fill(visited_vertex.begin(), visited_vertex.end(), false);
            std::fill(visited_node.begin(), visited_node.end(), false);
        }
    }
}

void GenerateQueryPlan::traversal_node(const Graph *query_graph, Edges***edge_matrix,
                                       std::vector<nd_tree_node> &density_tree,
                                       std::vector<uint32_t> &vertex_order, std::vector<uint32_t> &node_order,
                                       vector<bool> &visited_vertex, std::vector<bool> &visited_node,
                                       nd_tree_node &cur_node, std::unordered_set<uint32_t> &extendable_vertex) {
    node_order.push_back(cur_node.id_);
    visited_node[cur_node.id_] = true;

    // Check if all vertices has been visited.
    bool all_vertices_visited = true;
    for (auto u : cur_node.vertices_) {
        if (!visited_vertex[u]) {
            all_vertices_visited = false;
            break;
        }
    }

    if (all_vertices_visited) {
        if (cur_node.parent_ != -1 && !visited_node[cur_node.parent_])
            traversal_node(query_graph, edge_matrix, density_tree, vertex_order, node_order, visited_vertex, visited_node,
                           density_tree[cur_node.parent_], extendable_vertex);
        return;
    }

    std::vector<uint32_t> prev(query_graph->getVerticesCount());
    std::vector<double> dist(query_graph->getVerticesCount());

    // True means that there exists a child that is not visited.
    bool flag = true;
    while (flag) {
        double max_connectivity = 1.5;
        uint32_t selected_child = 0;
        flag = false;

        for (auto child_id : cur_node.children_) {
            if (!visited_node[child_id]) {
                double connectivity = connectivity_common_neighbors(visited_vertex, density_tree[child_id]);
                if (connectivity > max_connectivity
                    || (connectivity == max_connectivity &&
                        density_tree[child_id].density_ > density_tree[selected_child].density_)
                    || (connectivity == max_connectivity &&
                        density_tree[child_id].density_ > density_tree[selected_child].density_ &&
                        density_tree[child_id].vertices_.size() > density_tree[selected_child].vertices_.size())) {
                    max_connectivity = connectivity;
                    selected_child = child_id;
                }

                flag = true;
            }
        }

        if (flag) {
            if (max_connectivity < 2) {
                connectivity_shortest_path(query_graph, edge_matrix, vertex_order, visited_vertex, cur_node, prev, dist);
                max_connectivity = 0;
                uint32_t selected_vertex = 0;
                for (auto child_id : cur_node.children_) {
                    if (!visited_node[child_id]) {
                        double local_min_dist = std::numeric_limits<double>::max();
                        uint32_t local_selected_vertex = 0;
                        for (auto u : density_tree[child_id].vertices_) {
                            if (local_min_dist > dist[u]) {
                                local_min_dist = dist[u];
                                local_selected_vertex = u;
                            }
                        }

                        double connectivity = 1.0 / (local_min_dist + 1.0);
                        if (connectivity > max_connectivity) {
                            max_connectivity = connectivity;
                            selected_child = child_id;
                            selected_vertex = local_selected_vertex;
                        }
                    }
                }

                // Construct the path.
                std::vector<uint32_t> vertex_in_shortest_path;
                while (prev[selected_vertex] != selected_vertex) {
                    vertex_in_shortest_path.push_back(selected_vertex);
                    selected_vertex = prev[selected_vertex];
                }

                for (auto rb = vertex_in_shortest_path.rbegin(); rb != vertex_in_shortest_path.rend(); rb++) {
                    if (!visited_vertex[*rb]) {
                        vertex_order.push_back(*rb);
                        visited_vertex[*rb] = true;
                        update_extendable_vertex(query_graph, *rb, extendable_vertex, visited_vertex);

                        {
                            uint32_t bn_cnt_threshold = 2;
                            if (density_tree[selected_child].r_ > bn_cnt_threshold) {
                                bn_cnt_threshold = density_tree[selected_child].r_;
                            }

                            greedy_expand(query_graph, edge_matrix, vertex_order, visited_vertex, extendable_vertex,
                                          bn_cnt_threshold);
                        }
                    }
                }
            }

            traversal_node(query_graph, edge_matrix, density_tree, vertex_order, node_order, visited_vertex, visited_node,
                           density_tree[selected_child], extendable_vertex);

        }
    }


    if (vertex_order.empty()) {
        uint32_t min_relation_size = std::numeric_limits<uint32_t>::max();
        uint32_t max_degree_sum = 0;
        uint32_t first_vertex = 0;
        uint32_t second_vertex = 0;
        for (auto u : cur_node.vertices_) {
            uint32_t u_deg = query_graph->getVertexDegree(u);
            for (auto v : cur_node.vertices_) {
                uint32_t v_deg = query_graph->getVertexDegree(v);
                uint32_t local_degree_sum = u_deg + v_deg;
                if (u < v && query_graph->checkEdgeExistence(u, v)) {
                    if (local_degree_sum > max_degree_sum ||
                            (local_degree_sum == max_degree_sum &&
                                // storage->get_edge_relation_cardinality(u, v) < min_relation_size)) {
                                edge_matrix[u][v]->edge_count_ < min_relation_size)) {
                        max_degree_sum = local_degree_sum;
                        min_relation_size = edge_matrix[u][v]->edge_count_;
                        first_vertex = u;
                        second_vertex = v;

                        if (u_deg < v_deg) {
                            swap(first_vertex, second_vertex);
                        }
                    }
                }
            }
        }

        vertex_order.push_back(first_vertex);
        visited_vertex[first_vertex] = true;
        update_extendable_vertex(query_graph, first_vertex, extendable_vertex, visited_vertex);

        vertex_order.push_back(second_vertex);
        visited_vertex[second_vertex] = true;
        update_extendable_vertex(query_graph, second_vertex, extendable_vertex, visited_vertex);
    }

    while (!all_vertices_visited) {
        all_vertices_visited = true;

        uint32_t max_bn = 0;
        uint32_t max_degree = 0;
        uint32_t min_candidate_vertex_num = std::numeric_limits<uint32_t>::max();
        uint32_t selected_vertex = 0;
        for (auto u : cur_node.vertices_) {
            if (!visited_vertex[u]) {
                uint32_t nbrs_cnt;
                const uint32_t *nbrs = query_graph->getVertexNeighbors(u, nbrs_cnt);
                uint32_t local_bn_cnt = 0;
                uint32_t u_deg = query_graph->getVertexDegree(u);
                // uint32_t local_candidate_vertex_num = storage->num_candidates_[u];
                uint32_t local_candidate_vertex_num = edge_matrix[u][nbrs[0]]->vertex_count_;

                for (uint32_t i = 0; i < nbrs_cnt; ++i) {
                    uint32_t v = nbrs[i];

                    if (visited_vertex[v]) {
                        local_bn_cnt += 1;
                    }
                }

                if (local_bn_cnt > max_bn ||
                        (local_bn_cnt == max_bn && local_candidate_vertex_num < min_candidate_vertex_num) ||
                    (local_bn_cnt == max_bn &&local_candidate_vertex_num == min_candidate_vertex_num && u_deg > max_degree)) {
                    selected_vertex = u;
                    min_candidate_vertex_num = local_candidate_vertex_num;
                    max_bn = local_bn_cnt;
                    max_degree = u_deg;
                }

                all_vertices_visited = false;
            }
        }

        if (!all_vertices_visited) {
            vertex_order.push_back(selected_vertex);
            visited_vertex[selected_vertex] = true;
            update_extendable_vertex(query_graph, selected_vertex, extendable_vertex, visited_vertex);
        }
    }
    {
        uint32_t bn_cnt_threshold = 2;
        if (cur_node.parent_ != -1 && density_tree[cur_node.parent_].r_ > bn_cnt_threshold) {
            bn_cnt_threshold = density_tree[cur_node.parent_].r_;
        }

        greedy_expand(query_graph, edge_matrix, vertex_order, visited_vertex, extendable_vertex, bn_cnt_threshold);
    }
    if (cur_node.parent_ != -1 && !visited_node[cur_node.parent_])
        traversal_node(query_graph, edge_matrix, density_tree, vertex_order, node_order, visited_vertex, visited_node,
                       density_tree[cur_node.parent_], extendable_vertex);
}

double GenerateQueryPlan::connectivity_common_neighbors(std::vector<bool> &visited_vertex, nd_tree_node &child_node) {
    double connectivity = 1;
    for (auto u : child_node.vertices_) {
        if (visited_vertex[u]) {
            connectivity += 1;
        }
    }
    return connectivity;
}

void GenerateQueryPlan::connectivity_shortest_path(const Graph *query_graph, Edges***edge_matrix,
                                                   std::vector<uint32_t> &vertex_order,
                                                   std::vector<bool> &visited_vertex, nd_tree_node &cur_node,
                                                   std::vector<uint32_t> &prev, std::vector<double> &dist) {
    std::fill(prev.begin(), prev.end(), 0);
    std::fill(dist.begin(), dist.end(), std::numeric_limits<double>::max());
    std::vector<bool> in_cur_node(query_graph->getVerticesCount(), false);

    priority_queue<std::pair<double, uint32_t>, std::vector<std::pair<double, uint32_t>>, std::greater<std::pair<double, uint32_t>>> min_pq;

    for (auto u : vertex_order) {
        if (visited_vertex[u]) {
            prev[u] = u;
        }
        min_pq.push(std::make_pair(1, u));
    }

    for (auto u : cur_node.vertices_) {
        in_cur_node[u] = true;
    }

    while (!min_pq.empty()) {
        auto top_element = min_pq.top();
        min_pq.pop();

        uint32_t u = top_element.second;
        double dist_in_pq = top_element.first;

        // If the distance is updated, then skip the element.
        if (dist_in_pq > dist[u])
            continue;

        uint32_t nbrs_cnt;
        const uint32_t* nbrs = query_graph->getVertexNeighbors(u, nbrs_cnt);

        for (uint32_t i = 0; i < nbrs_cnt; ++i) {
            uint32_t v = nbrs[i];
            if (in_cur_node[v] && !visited_vertex[v]) {
                // double distance = storage->num_candidates_[v] > 3 ? storage->num_candidates_[v] : 3;
                double distance = edge_matrix[v][u]->vertex_count_ > 3 ? edge_matrix[v][u]->vertex_count_ : 3;
                double updated_dist = dist_in_pq * log(distance);

                if (updated_dist < dist[v]) {
                    dist[v] = updated_dist;
                    prev[v] = u;
                    min_pq.push(std::make_pair(updated_dist, v));
                }
            }
        }
    }
}

void GenerateQueryPlan::greedy_expand(const Graph *query_graph, Edges***edge_matrix, std::vector<uint32_t> &vertex_order,
                                      std::vector<bool> &visited_vertex,
                                      std::unordered_set<uint32_t> &extendable_vertex,
                                      uint32_t bn_cnt_threshold) {
    bool updated = true;

    while (updated) {
        updated = false;

        uint32_t max_bn_cnt = 0;
        uint32_t max_degree = 0;
        uint32_t min_candidate_vertex_num = std::numeric_limits<uint32_t>::max();
        uint32_t selected_vertex = 0;

        for (auto u : extendable_vertex) {
            uint32_t nbrs_cnt;
            const uint32_t* nbrs = query_graph->getVertexNeighbors(u, nbrs_cnt);

            uint32_t bn_cnt = 0;
            // uint32_t local_candidate_vertex_num = storage->num_candidates_[u];
            uint32_t local_candidate_vertex_num = edge_matrix[u][nbrs[0]]->vertex_count_;
            uint32_t u_deg = query_graph->getVertexDegree(u);

            for (uint32_t i = 0; i < nbrs_cnt; ++i) {
                uint32_t v = nbrs[i];
                if (visited_vertex[v]) {
                    bn_cnt += 1;
                }
            }

            if (bn_cnt > max_bn_cnt || (bn_cnt == max_bn_cnt && local_candidate_vertex_num < min_candidate_vertex_num)
                || (bn_cnt == max_bn_cnt && local_candidate_vertex_num == min_candidate_vertex_num && u_deg > max_degree)) {
                max_bn_cnt = bn_cnt;
                max_degree = u_deg;
                min_candidate_vertex_num = local_candidate_vertex_num;
                selected_vertex = u;
            }
        }

        if (max_bn_cnt >= bn_cnt_threshold) {
            vertex_order.push_back(selected_vertex);
            visited_vertex[selected_vertex] = true;
            update_extendable_vertex(query_graph, selected_vertex, extendable_vertex, visited_vertex);
            updated = true;
        }
    }
}

void GenerateQueryPlan::update_extendable_vertex(const Graph *query_graph, uint32_t u,
                                                 std::unordered_set<uint32_t> &extendable_vertex,
                                                 vector<bool> &visited_vertex) {
    uint32_t nbrs_cnt;
    const uint32_t* nbrs = query_graph->getVertexNeighbors(u, nbrs_cnt);

    for (uint32_t i = 0; i < nbrs_cnt; ++i) {
        uint32_t v = nbrs[i];
        if (!visited_vertex[v]) {
            extendable_vertex.insert(v);
        }
    }

    if (extendable_vertex.count(u) != 0) {
        extendable_vertex.erase(u);
    }
}

void GenerateQueryPlan::select_vertex_order(const Graph *query_graph, std::vector<std::vector<uint32_t>> &vertex_orders,
                                            std::vector<std::vector<uint32_t>> &node_orders) {

    std::vector<std::vector<uint32_t>> selected_vertex_order;
    uint32_t max_utility_value = 0;
    uint32_t selected_vertex_order_id = 0;
    for (uint32_t i = 0; i < vertex_orders.size(); ++i) {
        std::vector<uint32_t> bn_cnt_list;
        uint32_t utility_value = query_plan_utility_value(query_graph, vertex_orders[i], bn_cnt_list);

        if (utility_value > max_utility_value) {
            max_utility_value = utility_value;
            selected_vertex_order_id = i;
        }
    }

    selected_vertex_order.emplace_back(vertex_orders[selected_vertex_order_id]);
    vertex_orders.swap(selected_vertex_order);
}

uint32_t GenerateQueryPlan::query_plan_utility_value(const Graph *query_graph, std::vector<uint32_t> &vertex_order,
                                                     std::vector<uint32_t> &bn_cn_list) {
    uint32_t n = query_graph->getVerticesCount();
    std::vector<bool> visited(n, false);

    for (auto u : vertex_order) {
        visited[u] = true;
        uint32_t nbrs_cnt;
        const uint32_t *nbrs = query_graph->getVertexNeighbors(u, nbrs_cnt);

        uint32_t bn_cnt = 0;
        for (uint32_t i = 0; i < nbrs_cnt; ++i) {
            if (visited[nbrs[i]]) {
                bn_cnt += 1;
            }
        }

        bn_cn_list.push_back(bn_cnt);
    }

    uint32_t sum = 0;

    for (uint32_t i = 0; i < n; ++i) {
        for (uint32_t j = 0; j <= i; ++j) {
            sum += bn_cn_list[j];
        }
    }

    return sum;
}
