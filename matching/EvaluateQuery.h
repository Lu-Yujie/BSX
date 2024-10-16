#ifndef SUBGRAPHMATCHING_EVALUATEQUERY_H
#define SUBGRAPHMATCHING_EVALUATEQUERY_H

#include "utility/bsx/bsx.h"
#include <vector>
#include <queue>
#include <unordered_set>
#include <bitset>
#include <gmp.h>

class EvaluateQuery {
public:
    static void
    BS1Engine(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                  ui *candidates_count, ui *order, ui *pivot, size_t output_limit_num, size_t &call_count, mpz_t embedding_cnt);

    static void
    BSXEngine(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
               ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, mpz_t embedding_cnt);

private:
    static void generateBN(const Graph *query_graph, ui *order, ui *pivot, ui **&bn, ui *&bn_count);
    static void allocateBuffer(const Graph *query_graph, const Graph *data_graph, ui *candidates_count, ui *&idx,
                                   ui *&idx_count, ui *&embedding, ui *&idx_embedding, ui *&temp_buffer,
                                   ui **&valid_candidate_idx, bool *&visited_vertices);
    static void releaseBuffer(ui q_num, ui *idx, ui *idx_count, ui *embedding, ui *idx_embedding,
                                  ui *temp_buffer, ui **valid_candidate_idx, bool *visited_vertices, ui **bn, ui *bn_count);

    static void generateValidCandidateIndex(const Graph *data_graph, ui depth, ui *embedding, ui *idx_embedding,
                                            ui *idx_count, ui **valid_candidate_index, Edges ***edge_matrix,
                                            bool *visited_vertices, ui **bn, ui *bn_cnt, ui *order, ui *pivot,
                                            ui **candidates, const Graph *query_graph);

    static void bsxMaxCoverOrder(const Graph *graph, ui*& order, ui& num_cover, ui *candidates_count);

    static void bsxDeRefine(BSXIndex& index);

    static VertexID bsxGenNxtU(BSXIndex& index, VertexID* order, ui depth, ui num_cover);

    static bool bsxCheckTermination(ui num, VertexID* indep, std::stack<ui>*valid_cnt);

    static bool bsxGenIndepValidCans(ui indep_num, const VertexID* indep, BSXIndex& index, std::vector<std::vector<VertexID>>& cans);

    static ui bsxRefine(BSXIndex& index, VertexID u);

    static void bsxComEqBatch(BSXIndex& index, VertexID u);

    static ui bsxSepDiff(std::vector<VertexID> &v_cans, const ui *indep_con_cnt, int forward_idx, int backward_idx);

    static void bsxEnumerate4Parts(ui **&sep_flags, const VertexID* nodes, ui num_nodes, std::vector<std::vector<VertexID>>& cans, bool *&visited_v, mpz_t cur_cnt);

    static void bsxComEqBatchDirect(BSXIndex& index, VertexID u, std::vector<ui>& idxs);

    static void bsxGenResult(ui indep_num, const VertexID* indep, BSXIndex& index);

};


#endif //SUBGRAPHMATCHING_EVALUATEQUERY_H
