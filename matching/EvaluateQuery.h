#ifndef SUBGRAPHMATCHING_EVALUATEQUERY_H
#define SUBGRAPHMATCHING_EVALUATEQUERY_H

#include "QFilter.h"
#include "rapidMatch/primitive/projection.h"
#include "rapidMatch/relation/catalog.h"
#include "bsx/bsx.h"
#include "FiPE/FiPE.h"
#include <vector>
#include <queue>
#include <bitset>
#include <unordered_set>

class EvaluateQuery {
public:
    static bool BS1Engine(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                               ui *candidates_count, ui *order, ui *pivot, uint64_t output_limit_num, uint64_t &call_cnt,
                               mpz_t embedding_cnt, int64_t& time_limit);

    static bool RMEngine(const Graph *query_graph, const Graph *data_graph, catalog*&storage, Edges***edge_matrix,
                               ui **candidates, ui *candidates_count, ui *order, uint64_t output_limit_num, uint64_t &call_cnt,
                               mpz_t embedding_cnt, int64_t& time_limit);

    static bool KSSEngine(const Graph *query_graph, const Graph *data_graph, Edges***edge_matrix, ui **candidates, ui *candidates_count,
                                 ui *order, uint64_t output_limit_num, uint64_t &call_cnt, mpz_t embedding_cnt, int64_t& time_limit);

    static bool BSXEngine(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                           uint64_t output_limit_num, uint64_t &call_cnt, mpz_t embedding_cnt, int64_t& time_limit);

    static bool FiPEEngine(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                              ui *candidates_count, size_t output_limit_num, size_t &call_count, mpz_t embedding_cnt, int64_t& time_limit);

private:
    static void generateBN(const Graph *query_graph, ui *order, ui *pivot, ui **&bn, ui *&bn_count);

    static void allocateBuffer(const Graph *query_graph, const Graph *data_graph, ui *candidates_count, ui *&idx,
                                   ui *&idx_count, ui *&embedding, ui *&idx_embedding, ui *&temp_buffer,
                                   ui **&valid_candidate_idx, bool *&visited_vertices);

    static void releaseBuffer(ui q_num, ui *idx, ui *idx_count, ui *embedding, ui *idx_embedding,
                                  ui *temp_buffer, ui **valid_candidate_idx, bool *visited_vertices, ui **bn, ui *bn_count);

    static void exploreGenValidCanIdx(const Graph *data_graph, ui depth, ui *embedding, ui *idx_embedding,
                                            ui *idx_count, ui **valid_candidate_index, Edges ***edge_matrix,
                                            bool *visited_vertices, ui **bn, ui *bn_cnt, ui *order, ui *pivot,
                                            ui **candidates, const Graph *query_graph);

    static void convertCans2Catalog(const Graph *query_graph, ui **candidates, Edges ***edge_matrix, catalog *storage);

    static void convert_to_encoded_relation(catalog *storage, ui *order);

    static void convert_to_encoded_relation(catalog *storage, uint32_t u, uint32_t v);

    static void convert_to_hash_relation(catalog *storage, uint32_t u, uint32_t v);

    static void convert_encoded_relation_to_sparse_bitmap(catalog *storage, ui*order);

    static void updateShell2Kernel(const Graph *query_graph, VertexID u, ui* shell2kernel, bool* kos, std::vector<VertexID> & update);

    static void restoreShell2Kernel(const Graph *query_graph, VertexID u, ui* shell2kernel, bool* kos);

    static void kssComValidCans(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count,
                      ui**valid_cans, ui*valid_cans_count, ui* embedding, VertexID u, bool* visited_u, bool * visited_v);

    static bool kssGenResult(ui shell_num, ui* shell, ui** valid_cans, ui* valid_cans_count, bool * visited_v,
                                         mpz_t embedding_cnt, int64_t& time_limit);

    static bool kssGenResultImpl(ui depth, ui shell_num, ui* shell, ui** valid_cans, ui* valid_cans_count, bool * visited_v,
                                              mpz_t embedding_cnt, int64_t& time_limit);

    static void maxCoverOrder(const Graph *graph, ui*& order, ui& num_cover, ui *candidates_count);

    static void bsxDeRefine(BSXIndex& index);

    static VertexID bsxGenNxtU(BSXIndex& index, VertexID* order, ui depth, ui num_cover);

    static bool bsxCheckTermination(ui num, VertexID* indep, std::stack<ui>*valid_cnt);

    static bool bsxGenIndepValidCans(ui indep_num, const VertexID* indep, BSXIndex& index, std::vector<std::vector<VertexID>>& cans);

    static ui bsxRefine(BSXIndex& index, VertexID u);

    static void bsxComEqBatch(BSXIndex& index, VertexID u);

    static ui sepDiff(std::vector<VertexID> &v_cans, const ui *indep_con_cnt, int forward_idx, int backward_idx);

    static void enum4Parts(ui **&sep_flags, const VertexID* nodes, ui num_nodes, std::vector<std::vector<VertexID>>& cans, bool *&visited_v, mpz_t cur_cnt);

    static void bsxComEqBatchDirect(BSXIndex& index, VertexID u, std::vector<ui>& idxs);

    static void bsxGenResult(ui indep_num, const VertexID* indep, BSXIndex& index);

    static void FiPEEnum(FiPEIndex& index);

    static bool comSharedDis(FiPEIndex& index, ui depth);

    static bool comSharedCon(FiPEIndex& index, ui depth);

    static bool comCurSpace(FiPEIndex& index, ui depth);

    static bool setCurSpaceCon(FiPEIndex& index, ui depth);

    static void clearCurSpaceCon(FiPEIndex& index, ui depth);

    static void setCurSpaceDis(FiPEIndex& index, ui depth);

    static void clearCurSpaceDis(FiPEIndex& index, ui depth);

    static void comSub(FiPEIndex& index, ui depth);

    static void comEdgeSubCon(FiPEIndex& index, ui depth);

    static void comEdgeSubDis(FiPEIndex& index, ui depth);

    static bool comStartCans(FiPEIndex& index);

    inline static void splitCans(FiPEIndex& index, ui depth);
    inline static bool nxtSubCans(FiPEIndex& index, ui depth);

};


#endif //SUBGRAPHMATCHING_EVALUATEQUERY_H
