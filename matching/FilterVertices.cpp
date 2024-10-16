#include "FilterVertices.h"
#include "GenerateFilteringPlan.h"
#include <memory.h>
#include <utility/graphoperations.h>
#include <vector>
#include <algorithm>

/**
 * label and degree filtering
 * has nothing to do with edge
*/
bool
FilterVertices::LDFFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count) {
    allocateBuffer(data_graph, query_graph, candidates, candidates_count);

    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        LabelID label = query_graph->getVertexLabel(i);
        ui degree = query_graph->getVertexDegree(i);

        ui data_vertex_num;
        const ui* data_vertices = data_graph->getVerticesByLabel(label, data_vertex_num);

        for (ui j = 0; j < data_vertex_num; ++j) {
            ui data_vertex = data_vertices[j];
            if (data_graph->getVertexDegree(data_vertex) >= degree) {
                candidates[i][candidates_count[i]++] = data_vertex;
            }
        }

        if (candidates_count[i] == 0) {
            return false;
        }
    }

    return true;
}

/**
 * neighbor label frequency filtering
 * with edge label: add edge label checking when check neighbor labels
*/
bool
FilterVertices::NLFFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count) {
    allocateBuffer(data_graph, query_graph, candidates, candidates_count);

    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        VertexID query_vertex = i;
        computeCandidateWithNLF(data_graph, query_graph, query_vertex, candidates_count[query_vertex], candidates[query_vertex]);

        if (candidates_count[query_vertex] == 0) {
            return false;
        }
    }

    return true;
}

// 
bool
FilterVertices::CFLFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count,
                          ui *&order, TreeNode *&tree) {
    allocateBuffer(data_graph, query_graph, candidates, candidates_count);
    int level_count;
    ui* level_offset;
    GenerateFilteringPlan::generateCFLFilterPlan(data_graph, query_graph, tree, order, level_count, level_offset);

    VertexID start_vertex = order[0];
    computeCandidateWithNLF(data_graph, query_graph, start_vertex, candidates_count[start_vertex], candidates[start_vertex]);

    ui* updated_flag = new ui[data_graph->getVerticesCount()];
    ui* flag = new ui[data_graph->getVerticesCount()];
    std::fill(flag, flag + data_graph->getVerticesCount(), 0);

    // Top-down generation.
    for (int i = 1; i < level_count; ++i) {
        // Forward generation.
        for (int j = level_offset[i]; j < level_offset[i + 1]; ++j) {
            VertexID query_vertex = order[j];
            TreeNode& node = tree[query_vertex];
            generateCandidates(data_graph, query_graph, query_vertex, node.bn_, node.bn_count_, candidates, candidates_count, flag, updated_flag);
        }

        // Backward prune.
        for (int j = level_offset[i + 1] - 1; j >= level_offset[i]; --j) {
            VertexID query_vertex = order[j];
            TreeNode& node = tree[query_vertex];

            if (node.fn_count_ > 0) {
                pruneCandidates(data_graph, query_graph, query_vertex, node.fn_, node.fn_count_, candidates, candidates_count, flag, updated_flag);
            }
        }
    }

    // Bottom-up refinement.
    for (int i = level_count - 2; i >= 0; --i) {
        for (int j = level_offset[i]; j < level_offset[i + 1]; ++j) {
            VertexID query_vertex = order[j];
            TreeNode& node = tree[query_vertex];

            if (node.under_level_count_ > 0) {
                pruneCandidates(data_graph, query_graph, query_vertex, node.under_level_, node.under_level_count_, candidates, candidates_count, flag, updated_flag);
            }
        }
    }


    compactCandidates(candidates, candidates_count, query_graph->getVerticesCount());

    delete[] updated_flag;
    delete[] flag;
    return isCandidateSetValid(candidates, candidates_count, query_graph->getVerticesCount());
}

bool
FilterVertices::DPisoFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count,
                            ui *&order, TreeNode *&tree) {
    if (!LDFFilter(data_graph, query_graph, candidates, candidates_count))
        return false;

    GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, tree, order);

    ui q_num = query_graph->getVerticesCount();
    ui* updated_flag = new ui[data_graph->getVerticesCount()];
    ui* flag = new ui[data_graph->getVerticesCount()];
    std::fill(flag, flag + data_graph->getVerticesCount(), 0);

    // The number of refinement is k. According to the original paper, we set k as 3.
    for (ui k = 0; k < 3; ++k) {
        if (k % 2 == 0) {
            for (int i = 1; i < q_num; ++i) {
                VertexID query_vertex = order[i];
                TreeNode& node = tree[query_vertex];
                pruneCandidates(data_graph, query_graph, query_vertex, node.bn_, node.bn_count_, candidates, candidates_count, flag, updated_flag);
            }
        }
        else {
            for (int i = q_num - 2; i >= 0; --i) {
                VertexID query_vertex = order[i];
                TreeNode& node = tree[query_vertex];
                pruneCandidates(data_graph, query_graph, query_vertex, node.fn_, node.fn_count_, candidates, candidates_count, flag, updated_flag);
            }
        }
    }

    compactCandidates(candidates, candidates_count, query_graph->getVerticesCount());

    delete[] updated_flag;
    delete[] flag;
    return isCandidateSetValid(candidates, candidates_count, query_graph->getVerticesCount());
}
void FilterVertices::allocateBuffer(const Graph *data_graph, const Graph *query_graph, ui **&candidates,
                                    ui *&candidates_count) {
    ui q_num = query_graph->getVerticesCount();
    ui candidates_max_num = data_graph->getGraphMaxLabelFrequency();

    candidates_count = new ui[q_num];
    memset(candidates_count, 0, sizeof(ui) * q_num);

    candidates = new ui*[q_num];

    for (ui i = 0; i < q_num; ++i) {
        candidates[i] = new ui[candidates_max_num];
    }
}

void FilterVertices::compactCandidates(ui **&candidates, ui *&candidates_count, ui q_num) {
    for (ui i = 0; i < q_num; ++i) {
        VertexID query_vertex = i;
        ui next_position = 0;
        for (ui j = 0; j < candidates_count[query_vertex]; ++j) {
            VertexID data_vertex = candidates[query_vertex][j];

            if (data_vertex != (ui)-1) {
                candidates[query_vertex][next_position++] = data_vertex;
            }
        }

        candidates_count[query_vertex] = next_position;
    }
}

bool FilterVertices::isCandidateSetValid(ui **&candidates, ui *&candidates_count, ui q_num) {
    for (ui i = 0; i < q_num; ++i) {
        if (candidates_count[i] == 0)
            return false;
    }
    return true;
}

void
FilterVertices::computeCandidateWithNLF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                               ui &count, ui *buffer) {
    LabelID label = query_graph->getVertexLabel(query_vertex);
    ui degree = query_graph->getVertexDegree(query_vertex);
#if OPTIMIZED_VLABELED_GRAPH == 1
    auto u_nlf = query_graph->getVertexNLF(query_vertex);
#endif
    ui data_vertex_num;
    const ui* data_vertices = data_graph->getVerticesByLabel(label, data_vertex_num);
    count = 0;
    for (ui j = 0; j < data_vertex_num; ++j) {
        ui data_vertex = data_vertices[j];
        if (data_graph->getVertexDegree(data_vertex) >= degree) {

            // NLF check
#if OPTIMIZED_VLABELED_GRAPH == 1
            auto v_nlf = data_graph->getVertexNLF(data_vertex);
            if (v_nlf->size() >= u_nlf->size()) {
                bool is_valid = true;
                for (auto element : *u_nlf) {
                    auto iter = v_nlf->find(element.first);
                    if (iter == v_nlf->end() || iter->second < element.second) {
                        is_valid = false;
                        break;
                    }
                }

                if (is_valid) {
                    if (buffer != nullptr) {
                        buffer[count] = data_vertex;
                    }
                    count += 1;
                }
            }
#else
            if (buffer != nullptr) {
                buffer[count] = data_vertex;
            }
            count += 1;
#endif
        }
    }

}

void FilterVertices::computeCandidateWithLDF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                             ui &count, ui *buffer) {
    LabelID label = query_graph->getVertexLabel(query_vertex);
    ui degree = query_graph->getVertexDegree(query_vertex);
    count = 0;
    ui data_vertex_num;
    const ui* data_vertices = data_graph->getVerticesByLabel(label, data_vertex_num);

    if (buffer == nullptr) {
        for (ui i = 0; i < data_vertex_num; ++i) {
            VertexID v = data_vertices[i];
            if (data_graph->getVertexDegree(v) >= degree) {
                count += 1;
            }
        }
    }
    else {
        for (ui i = 0; i < data_vertex_num; ++i) {
            VertexID v = data_vertices[i];
            if (data_graph->getVertexDegree(v) >= degree) {
                buffer[count++] = v;
            }
        }
    }
}

void FilterVertices::generateCandidates(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                       VertexID *pivot_vertices, ui pivot_vertices_count, VertexID **candidates,
                                       ui *candidates_count, ui *flag, ui *updated_flag) {
    LabelID query_vertex_label = query_graph->getVertexLabel(query_vertex);
    ui query_vertex_degree = query_graph->getVertexDegree(query_vertex);
#if OPTIMIZED_VLABELED_GRAPH == 1
    auto u_nlf = query_graph->getVertexNLF(query_vertex);
#endif
    ui count = 0;
    ui updated_flag_count = 0;
    for (ui i = 0; i < pivot_vertices_count; ++i) {
        VertexID pivot_vertex = pivot_vertices[i];

        for (ui j = 0; j < candidates_count[pivot_vertex]; ++j) {
            VertexID v = candidates[pivot_vertex][j];

            if (v == (ui)-1)
                continue;
            ui v_nbrs_count;
            const VertexID* v_nbrs = data_graph->getVertexNeighbors(v, v_nbrs_count);

            for (ui k = 0; k < v_nbrs_count; ++k) {
                VertexID v_nbr = v_nbrs[k];
                LabelID v_nbr_label = data_graph->getVertexLabel(v_nbr);
                ui v_nbr_degree = data_graph->getVertexDegree(v_nbr);

                if (flag[v_nbr] == count && v_nbr_label == query_vertex_label && v_nbr_degree >= query_vertex_degree) {
                    flag[v_nbr] += 1;

                    if (count == 0) {
                        updated_flag[updated_flag_count++] = v_nbr;
                    }
                }
            }
        }

        count += 1;
    }

    for (ui i = 0; i < updated_flag_count; ++i) {
        VertexID v = updated_flag[i];
        if (flag[v] == count) {
            // NLF filter.
#if OPTIMIZED_VLABELED_GRAPH == 1
            auto v_nlf = data_graph->getVertexNLF(v);
            if (v_nlf->size() >= u_nlf->size()) {
                bool is_valid = true;
                for (auto element : *u_nlf) {
                    auto iter = v_nlf->find(element.first);
                    if (iter == v_nlf->end() || iter->second < element.second) {
                        is_valid = false;
                        break;
                    }
                }
                if (is_valid) {
                    candidates[query_vertex][candidates_count[query_vertex]++] = v;
                }
            }
#else
            candidates[query_vertex][candidates_count[query_vertex]++] = v;
#endif
        }
    }

    for (ui i = 0; i < updated_flag_count; ++i) {
        ui v = updated_flag[i];
        flag[v] = 0;
    }
}

void
FilterVertices::pruneCandidates(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                    VertexID *pivot_vertices, ui pivot_vertices_count, VertexID **candidates,
                                    ui *candidates_count, ui *flag, ui *updated_flag) {
    LabelID query_vertex_label = query_graph->getVertexLabel(query_vertex);
    ui query_vertex_degree = query_graph->getVertexDegree(query_vertex);

    ui count = 0;
    ui updated_flag_count = 0;
    for (ui i = 0; i < pivot_vertices_count; ++i) {
        VertexID pivot_vertex = pivot_vertices[i];

        for (ui j = 0; j < candidates_count[pivot_vertex]; ++j) {
            VertexID v = candidates[pivot_vertex][j];

            if (v == (ui)-1)
                continue;
            ui v_nbrs_count;
            auto v_nbrs = data_graph->getVertexNeighbors(v, v_nbrs_count);

            for (ui k = 0; k < v_nbrs_count; ++k) {
                VertexID v_nbr = v_nbrs[k];
                LabelID v_nbr_label = data_graph->getVertexLabel(v_nbr);
                ui v_nbr_degree = data_graph->getVertexDegree(v_nbr);

                if (flag[v_nbr] == count && v_nbr_label == query_vertex_label && v_nbr_degree >= query_vertex_degree) {
                    flag[v_nbr] += 1;

                    if (count == 0) {
                        updated_flag[updated_flag_count++] = v_nbr;
                    }
                }
            }
        }

        count += 1;
    }
    for (ui i = 0; i < candidates_count[query_vertex]; ++i) {
        ui v = candidates[query_vertex][i];
        if (v == (ui)-1)
            continue;

        if (flag[v] != count) {
            candidates[query_vertex][i] = (ui)-1;
        }
    }

    for (ui i = 0; i < updated_flag_count; ++i) {
        ui v = updated_flag[i];
        flag[v] = 0;
    }
}

void FilterVertices::sortCandidates(ui **candidates, ui *candidates_count, ui num) {
    for (ui i = 0; i < num; ++i) {
        std::sort(candidates[i], candidates[i] + candidates_count[i]);
    }
}
