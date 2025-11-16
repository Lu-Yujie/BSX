#ifndef LU_BSX_H
#define LU_BSX_H

/**
 * define data structure used in BSX query
*/
#include <stack>
#include <unordered_map>
#include <vector>
#include <bitset>
#include <gmp.h>
#include "graph/graph.h"
#include "pretty_print.h"
using namespace std;
typedef unsigned int ui;

class b_search{
public:
/**
 * universal binary search for ui*
*/
// small array: linear seach
static ui smallArraySearch(const ui* arr, ui size, ui target) {
    for (ui i = 0; i < size; ++i) {
        if (arr[i] == target)
            return i;
    }
    return (ui)-1;
}

// large array, use lower bound func
static ui largeArraySearch(const ui* arr, ui size, ui target) {
    auto ptr = lower_bound(arr, arr + size, target);
    if (ptr != arr + size && *ptr == target)
        return ptr - arr;
    else
        return (ui)-1;
}

static ui search(const ui* arr, ui size, ui target) {
    if (size <= 4) { // 可根据经验设定阈值
        return smallArraySearch(arr, size, target);
    }
    return largeArraySearch(arr, size, target);
}
};  // class b_search

struct SimInfo {
    bool enable_;
    ui** batches_idxs_;
    ui** sim_batches_;
    ui** offset_;
    ui* num_;
    ui** idxs_;  // indicate each eq_batch belongs to which sim_batch
    bool* valid_;
    ui* max_diff_;
    VertexID** diff_unbrs_;  // the u influences diff_unbrs, match to u, set sim_batch of diff_unbrs to invalid
    ui* diff_unbrs_cnt_;
    VertexID* judge_nbrs_;
    bool** same_or_diff_;  // diff->true, same&no-relation->false
    bool** batches_processed_;  // used to skip invalid batches
    ui q_num;
    SimInfo() { enable_ = false; }
    void init(const Graph* graph, ui max_candidates_num) {
        q_num = graph->getVerticesCount();
        batches_idxs_ = new ui*[q_num];
        sim_batches_ = new ui*[q_num];
        offset_ = new ui*[q_num];
        num_ = new ui[q_num];
        idxs_ = new ui*[q_num];
        valid_ = new bool[q_num];
        max_diff_ = new ui[q_num];
        diff_unbrs_ = new VertexID*[q_num];
        diff_unbrs_cnt_ = new ui[q_num];
        same_or_diff_ = new bool*[q_num];
        judge_nbrs_ = new VertexID[q_num];
        batches_processed_ = new bool*[q_num];
        memset(diff_unbrs_cnt_, 0, sizeof(ui)*q_num);
        memset(valid_, false, sizeof(bool)*q_num);
        for (ui i = 0; i < q_num; i++) {
            diff_unbrs_[i] = new VertexID[graph->getVertexDegree(i)];
            max_diff_[i] = (graph->getVertexDegree(i) - 1)*(1 - BSX_SIM_THRESHOLD * 0.01);
            if (max_diff_[i] == 0) continue;
            batches_idxs_[i] = new ui[max_candidates_num];
            sim_batches_[i] = new ui[max_candidates_num];
            offset_[i] = new ui[max_candidates_num + 1];
            idxs_[i] = new ui[max_candidates_num];
            batches_processed_[i] = new bool[max_candidates_num];
            same_or_diff_[i] = new bool[q_num];
        }
    }
    ~SimInfo() {
        if (!enable_) return;
        for (ui i = 0; i < q_num; i++) {
            delete[] diff_unbrs_[i];
            if (max_diff_[i] == 0) continue;
            delete[] batches_idxs_[i];
            delete[] sim_batches_[i];
            delete[] offset_[i];
            delete[] idxs_[i];
            delete[] batches_processed_[i];
            delete[] same_or_diff_[i];
        }
        delete[] batches_idxs_;
        delete[] sim_batches_;
        delete[] offset_;
        delete[] num_;
        delete[] idxs_;
        delete[] valid_;
        delete[] max_diff_;
        delete[] diff_unbrs_;
        delete[] diff_unbrs_cnt_;
        delete[] batches_processed_;
        delete[] same_or_diff_;
    }
};

/**structures used to store batch info
 * nodes: store batches which are seperated by offset
 * offset: indicate the position of each batch
 * cnt: incidate the number of each batch, 24-3-7
 *   1. add cnt because the some batches may become invalid along refinement, offset is not enough, 24-3-7
 *   2. stack after u is visited, 24-3-7
 *   3. valid_cans|valid_cnt is the same as cur_batch_nodes|cur_batch_cnt, if u is visited, 24-3-7
 *   4. relationship between valid_cans & batch_info
 *      4.1 both of them of u are the same when building batches of u
 *      4.2 valid_cans may shrink along matching, batch_info is solid once after grouping
 *      4.3 BSXRefine is based on valid_cans, ComBatch is based on batch_info
 *      4.4 ComBatch generates batch_info and then synchronizes to valid_cans
 *      4.5 for each next/backstracking step, push/pop batch_info&valid_cans of u, push/pop valid_cans of other influenced nodes
*/
struct BatchInfo {
    // cur_batch_nodes = nodes.top()+offset.top()[idx[depth]]
    // cur_batch_cnt = cnt.top()[idx[depth]]
    stack<VertexID*> nodes_;
    stack<ui*> offset_;  // start from 0
    stack<ui*> cnt_;     // # of each batch
    stack<ui> idx_;      // point to current batch
    stack<ui> num_;      // #batches
    stack<ui> maxCnt_;   // max # of each batches
    static SimInfo* sim_info;  // release memory at the end of engine

    BatchInfo() {}
    ~BatchInfo() {
        while(!nodes_.empty()) {
            delete[] nodes_.top();
            nodes_.pop();
            delete[] offset_.top();
            offset_.pop();
            delete[] cnt_.top();
            cnt_.pop();
        }
    }
    // current(idx'th) batch of u
    inline VertexID* cur_batch(ui& cnt) {
        // if (nodes_.empty()) return nullptr;
        cnt = cnt_.top()[idx_.top()];
        return nodes_.top()+offset_.top()[idx_.top()];
    }
    // stack new batch info
    inline void add() {
        nodes_.push(nullptr);
        offset_.push(nullptr);
        cnt_.push(nullptr);
        idx_.push(0);
        num_.push(0);
        maxCnt_.push(0);
    }
    // pop top batch info
    inline void pop() {
        delete[] nodes_.top();
        nodes_.pop();
        delete[] offset_.top();
        offset_.pop();
        delete[] cnt_.top();
        cnt_.pop();
        idx_.pop();
        num_.pop();
        maxCnt_.pop();
    }
    inline void print() {
        cout << "******************* batch info *******************" << endl;
        cout << "offset: ";
        for (ui i = 0; i < num_.top(); i++) {
            cout << offset_.top()[i] << ',';
        }
        cout << endl;
        cout << "cnt: ";
        for (ui i = 0; i < num_.top(); i++) {
            cout << cnt_.top()[i] << ',';
        }
        cout << endl;
        cout << "nodes: ";
        for (ui i = 0; i < num_.top(); i++) {
            for (ui j = 0; j < cnt_.top()[i]; j++) {
                cout << nodes_.top()[offset_.top()[i] + j] << ',';
            }
        }
        cout << "\n#batch: " << num_.top() << ", idx: " << idx_.top() << ", maxCnt: " << maxCnt_.top() << endl;
        cout << endl;
    }
};

/**
 * Embedding info
 * mapping between depth, u(query node), v(data node)
 * u2v: u->v, each u only match to one v
 * v2depth: v->depth, too much v, use map instead of array
 * depth2u: depth->u, use vector for dynamic tree height
*/
class Embedding{
public:
    VertexID* u2v;               // u->v, use the 1st v of batch
    vector<VertexID> depth2u;    // depth->u

    Embedding(ui cnt) {
        u2v = new VertexID[cnt];
        depth2u.reserve(cnt);
    }
    ~Embedding() {
        delete[] u2v;
    }
};

/**
 * new index structure, update in time
*/
class BSXIndex {
public:
    stack<Edges*>** index_;  // can be used to judge edge existence, index_[u_1][u_2].size() != 0
    stack<VertexID*>* valid_cans_;  // diff from batch!!!, used to locate can in index_ to get nbrs
    stack<ui>* valid_cnt_;  // #valid_cans_
    VertexID** index_cans_;  // used for index_, indicate valid_cans of index.top
    ui* index_cnt_;         //            ''             valid_cnt      ''
    const Graph* q_graph_;
    const Graph* d_graph_;
    // u.top().v -> (nbr_cnt, nbrs),  sorted
    // stack<unordered_map<VertexID, pair<ui, VertexID*>>*>* cached_uv_nbrs_;
    stack<vector<VertexID>> influenced_u_;  // influenced nodes in each layer, src_u is 1'th
    BatchInfo* batch_info;  // array[q_num], use uid as idx, because u may be grouped multi-times
    bool* visited_u;
    bool* visited_v;
    Embedding* embedding;
    ui num_cover_;
    ui* indep_con_cnt_;    // count the number of times each indep_cans may conflict
    ui** sep_flag_;        // seperate indep cans

    // temporary embeddings
    mpz_t level_embeddings_;  // #embeddings of one depth
    mpz_t label_embeddings_;  // #embeddings of one kind of label, for enumeration
    BSXIndex(const Graph*q_graph, const Graph*d_graph, Edges ***index, ui **cans, ui *cans_cnt, ui num_cover) {
        q_graph_ = q_graph;
        d_graph_ = d_graph;
        num_cover_ = num_cover;
        auto qnum = q_graph->getVerticesCount();
        auto dnum = d_graph->getVerticesCount();
        auto num_indep = qnum - num_cover;
        index_ = new stack<Edges*>*[qnum];
        batch_info = new BatchInfo[qnum];
        visited_u = new bool[qnum];
        memset(visited_u, false, sizeof(bool)* qnum);
        visited_v = new bool[dnum];
        memset(visited_v, false, sizeof(bool)* dnum);
        valid_cans_ = new stack<VertexID*>[qnum];
        valid_cnt_ = new stack<ui>[qnum];
        index_cans_ = new VertexID*[qnum];
        index_cnt_ = new VertexID[qnum];
        indep_con_cnt_ = new ui[dnum];
        memset(indep_con_cnt_, 0, sizeof(ui)*dnum);
        sep_flag_ = new ui*[qnum];
        embedding = new Embedding(qnum);
        for (ui i = 0; i < qnum; i++) {
            // each node has 3 seperation point, seperate into four parts
            // 1 upward conflict: seperate nodes based on whether it shows upward
            // 2 downward   ''  :    ''     ''     ''        ''       ''   downward
            // 4 parts: up-down,up-x,x-down,x-x
            sep_flag_[i] = new ui[3];
        }
        // cached_uv_nbrs_ = new stack<unordered_map<VertexID, pair<ui, VertexID*>>*>[qnum];
        for (ui i = 0; i < qnum; i++) {  // i->id of u
            index_[i] = new stack<Edges*>[qnum];
            for (ui j = 0; j < qnum; j++) {  // j->id of nbrs, no-nullptr->edge exists
                if (index[i][j] != nullptr) index_[i][j].push(index[i][j]);
            }
            valid_cans_[i].push(cans[i]);
            valid_cnt_[i].push(cans_cnt[i]);
            index_cans_[i] = cans[i];
            index_cnt_[i] = cans_cnt[i];
            // cached_uv_nbrs_[i].push(nullptr);
        }
        mpz_init(level_embeddings_);
        mpz_init(label_embeddings_);
    }

    ~BSXIndex() {
        auto qnum = q_graph_->getVerticesCount();
        auto num_indep = qnum - num_cover_;
        delete[] batch_info;
        delete[] visited_u;
        delete[] visited_v;
        delete[] valid_cnt_;
        delete[] index_cans_;
        delete[] index_cnt_;
        delete[] indep_con_cnt_;
        for (ui i = 0; i < num_indep; i++) {
            delete[] sep_flag_[i];
        }
        delete[] sep_flag_;
        delete embedding;
        for (ui i = 0; i < qnum; i++) {
            for (ui j = 0; j < qnum; j++) {  // j->id of nbrs, no-nullptr->edge exists
                while (index_[i][j].size() > 1) {
                    delete index_[i][j].top();
                    index_[i][j].pop();
                }
            }
            delete[] index_[i];
            while (valid_cans_[i].size() > 1) {
                delete[] valid_cans_[i].top();
                valid_cans_[i].pop();
            }
            valid_cans_[i].pop();
        }
        delete[] index_;
        delete[] valid_cans_;
        mpz_clear(level_embeddings_);
        mpz_clear(label_embeddings_);
    }

    // get Neighbors of v(can of u_1) from u_1 to u_2
    // used to generate valid cans, by union nbrs, 24-3-7
    const VertexID* getNeighbors(VertexID u_1, VertexID u_2, VertexID v, ui& nbrs_cnt) {
        auto v_idx = b_search::search(index_cans_[u_1], index_cnt_[u_1], v);
        if (v_idx == (ui)-1) {
            cout << "can't find " << v << " in index_cans_[" << u_1 << "]: ";
            for (ui i = 0; i < index_cnt_[u_1]; i++) cout << index_cans_[u_1][i] << ", ";
            cout << endl;
            exit(1);
        }
        auto& edges = *(index_[u_1][u_2].top());
        nbrs_cnt = edges.offset_[v_idx+1] - edges.offset_[v_idx];
        return edges.edge_ + edges.offset_[v_idx];
    }

    // get all Neighbors(valid) of can_v of u, sorted
    // used to seperate nodes, 24-3-7
    VertexID* getNeighbors(VertexID u, VertexID v, ui& nbrs_cnt) {
        // if (cached_uv_nbrs_[u].top() != nullptr) {
        //     auto nbrs_iter = cached_uv_nbrs_[u].top()->find(v);
        //     if (nbrs_iter != cached_uv_nbrs_[u].top()->end()) {
        //         nbrs_cnt = nbrs_iter->second.first;
        //         return nbrs_iter->second.second;
        //     }
        // }
        ui unbrs_cnt;
        set<VertexID> nbrs;
        auto unbrs = q_graph_->getVertexNeighbors(u, unbrs_cnt);
        for (ui i = 0; i < unbrs_cnt; i++ ) {
            auto& edges = *(index_[u][unbrs[i]].top());
            auto v_idx = b_search::search(valid_cans_[u].top(), valid_cnt_[u].top(), v);
            if (v_idx == (ui)-1) {
                exit(1);
            }
            for (ui k = edges.offset_[v_idx]; k < edges.offset_[v_idx+1]; k++) {
                nbrs.emplace(edges.edge_[k]);
            }
        }
        nbrs_cnt = nbrs.size();
        VertexID* nbrs_ptr = new ui[nbrs_cnt];
        ui nbr_idx = 0;
        for (auto& nbr : nbrs) nbrs_ptr[nbr_idx++] = nbr;
        // // cache the result
        // if (cached_uv_nbrs_[u].top() == nullptr) {
        //     cached_uv_nbrs_[u].top() = new unordered_map<VertexID, pair<ui, VertexID*>>;
        // }
        // cached_uv_nbrs_[u].top()->emplace(v, make_pair(nbrs_cnt, nbrs_ptr));
        return nbrs_ptr;
    }

    /**update the structure, check nothing
     * based on new valid_cans, valid_cnt & influneced
    */
    ui updateIndex(bool* influenced) {
        auto q_num = q_graph_->getVerticesCount();
        bool* nbr_updated = new bool[q_num];
        memset(nbr_updated, false, sizeof(bool)*q_num);
        ui* cans_idx = new ui[d_graph_->getVerticesCount()];
        memset(cans_idx, 0, sizeof(ui)*d_graph_->getVerticesCount());
        vector<VertexID> updated_cans_idx;  // used to recover cans_idx
        vector<ui> temp_edges;  // used to record & build edge.edge_

        for (ui u = 0; u < q_num; u++) {
            if (!influenced[u]) continue;
            nbr_updated[u] = true;

            // build cans_idx of u, and then build edges from nbr to u.
            ui updated_cans_idx_cnt = 0;
            for (ui i = 0; i < valid_cnt_[u].top(); i++) {
                VertexID v = valid_cans_[u].top()[i];
                cans_idx[v] = i + 1;
                updated_cans_idx.emplace_back(v);
                updated_cans_idx_cnt++;
            }

            // update all nbrs influenced
            ui nbrs_cnt;
            auto nbrs = q_graph_->getVertexNeighbors(u, nbrs_cnt);
            for (ui i = 0; i < nbrs_cnt; i++) {
                auto& nbr = nbrs[i];
                if (nbr_updated[nbr]) continue;

                // update edges between u & nbr
                auto nbr2u_edges = new Edges;
                nbr2u_edges->vertex_count_ = valid_cnt_[nbr].top();
                nbr2u_edges->offset_ = new ui[nbr2u_edges->vertex_count_ + 1];

                auto u2nbr_edges = new Edges;
                u2nbr_edges->vertex_count_ = valid_cnt_[u].top();
                u2nbr_edges->offset_ = new ui[u2nbr_edges->vertex_count_ + 1];
                fill(u2nbr_edges->offset_, u2nbr_edges->offset_ + u2nbr_edges->vertex_count_ + 1, 0);

                ui local_edge_count = 0;
                // build edges from nbr to u, then build edges from u to nbr(based on pre-result)
                for (ui j = 0; j < nbr2u_edges->vertex_count_; j++) {
                    VertexID v = valid_cans_[nbr].top()[j];
                    nbr2u_edges->offset_[j] = local_edge_count;
                    ui nbr_v_nbrs_cnt;
                    auto nbr_v_nbrs = getNeighbors(nbr, u, v, nbr_v_nbrs_cnt);  // error
                    for (ui k = 0; k < nbr_v_nbrs_cnt; k++) {
                        auto nbr_v_nbr = nbr_v_nbrs[k];
                        if (cans_idx[nbr_v_nbr] != 0) {
                            ui position = cans_idx[nbr_v_nbr];
                            u2nbr_edges->offset_[position] += 1;
                            temp_edges.emplace_back(nbr_v_nbr);
                            local_edge_count++;
                        }
                    }
                }
                nbr2u_edges->offset_[nbr2u_edges->vertex_count_] = local_edge_count;
                nbr2u_edges->edge_count_ = local_edge_count;
                nbr2u_edges->edge_ = new VertexID[local_edge_count];
                copy(temp_edges.begin(), temp_edges.end(), nbr2u_edges->edge_);
                temp_edges.clear();

                u2nbr_edges->edge_count_ = local_edge_count;
                u2nbr_edges->edge_ = new VertexID[local_edge_count];
                for (ui j = 1; j <= u2nbr_edges->vertex_count_; ++j) {
                    u2nbr_edges->offset_[j] += u2nbr_edges->offset_[j - 1];
                }

                // build edges from u2nbr, based on nbr2u, revert edge start&end
                for (ui j = 0; j < nbr2u_edges->vertex_count_; ++j) {
                    VertexID start = valid_cans_[nbr].top()[j];
                    for (ui k = nbr2u_edges->offset_[j]; k < nbr2u_edges->offset_[j + 1]; k++) {
                        VertexID end = nbr2u_edges->edge_[k];  // end is valid_cans of u
                        u2nbr_edges->edge_[u2nbr_edges->offset_[cans_idx[end] - 1]++] = start;
                    }
                }
                for (ui j = u2nbr_edges->vertex_count_; j >= 1; --j) {
                    u2nbr_edges->offset_[j] = u2nbr_edges->offset_[j - 1];
                }
                u2nbr_edges->offset_[0] = 0;

                // shrink edges, based on offset
                // can't shrink infact, because it may due to cascade reaction.
                //    if valid_cans change(indep node didn't propagate), then update its nbrs
                // if we just change u side, then cascade reaction won't happen
                //   This will work because the nbrs side will not be influenced if u is changed
                //   And benefited from the storage of VertexID in edges instead of v_idx, we
                //   do not need to update the edges.edge_ of nbrs&u
                // But, if valid_cans of u is changed, we need to re-compute its all nbrs
                //   as a final method, run some tests if there are enough time, 24-4-5
                //   TODO: write delete method and do tests

                // add new edges(nbrs info)
                index_[nbr][u].push(nbr2u_edges);
                index_[u][nbr].push(u2nbr_edges);
            }

            // recover cans_idx for next u
            for (ui i = 0; i < updated_cans_idx_cnt; i++) {
                cans_idx[updated_cans_idx[i]] = 0;
            }
            updated_cans_idx.clear();
        }
        delete[] nbr_updated;
        delete[] cans_idx;
        return 0;
    }

    /**
     * just update one side of index_(edges), allow hanging nodes
     * hanging nodes: which have been deleted in valid_cans in its nbrs
     *                because the node is invalid for the other nbrs of its nbrs
     * 
     * based on new valid_cans, valid_cnt & influneced
     * shrink valid_cans, if offset[i]-offset[i-1] == 0
     * return -1, if valid_cans[*] == 0
    */
    ui updateOneSide(bool* influenced) {
        cout << "not implemented func" << endl;
        exit(0);
        return 0;
    }
};

#endif
