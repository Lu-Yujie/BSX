#ifndef LU_FIPE_H
#define LU_FIPE_H

/**
 * define data structure used in FiPE query
 */
#include <stack>
#include <unordered_map>
#include <bitset>
#include <gmp.h>
#include "common.h"
#include "pretty_print.h"
#include "bsx/SetOp.h"
#include "timeOp.h"
#include "FiPE/IndepSet.h"
#include "FiPE/edgeEqu.h"
using namespace std;
typedef unsigned int ui;

/**
 * new index structure, update in time
 */
class FiPEIndex {
public:
  const Graph *q_graph_;
  const Graph *d_graph_;
  bool *visited_u;
  bool *visited_v;
  deque<vector<VertexID>> *valid_cans_;
  VertexID* order_;

  deque<Edges *> **index_; // can be used to judge edge existence, index_[u_1][u_2].size() != 0
  Embedding *embedding;
  FiPEIndep *indepInfo;
  SubInfo* subInfo_;
  SubCans* subCans_;  // u_id as idx

#ifdef ANALYZE_TIME
  static int64_t refine_time;
  static int64_t enumerate_time;
  static int64_t getNeighbors_time;
#endif

  FiPEIndex(const Graph *q_graph, const Graph *d_graph, Edges ***index, ui **cans, ui *cans_cnt,
               ui num_cover, VertexID* order) {
    q_graph_ = q_graph;
    d_graph_ = d_graph;
    auto qnum = q_graph->getVerticesCount();
    auto dnum = d_graph->getVerticesCount();
    visited_u = new bool[qnum];
    memset(visited_u, false, sizeof(bool) * qnum);
    visited_v = new bool[dnum];
    memset(visited_v, false, sizeof(bool) * dnum);
    order_ = order;

    valid_cans_ = new deque<vector<VertexID>>[qnum];
    index_ = new deque<Edges *> *[qnum];
    for (ui i = 0; i < qnum; i++) { // i->id of u
      index_[i] = new deque<Edges *>[qnum];
      for (ui j = 0; j < qnum; j++) { // j->id of nbrs, no-nullptr->edge exists
        if (index[i][j] != nullptr)
          index_[i][j].emplace_back(index[i][j]);
      }
      vector<VertexID> tmp_cans;
      tmp_cans.reserve(cans_cnt[i]);
      for (ui j = 0; j < cans_cnt[i]; j++)
        tmp_cans.emplace_back(cans[i][j]);
      valid_cans_[i].emplace_back(std::move(tmp_cans));
    }

    embedding = new Embedding(qnum);
    indepInfo = new FiPEIndep(qnum, dnum, num_cover, order);
    subInfo_ = new SubInfo[num_cover-1];
    for (ui i = 0; i < num_cover-1; i++) {
      visited_u[order[i]] = true;
      subInfo_[i].init(q_graph, order, i, visited_u, num_cover);
    }
    memset(visited_u, false, sizeof(bool) * qnum);
    subCans_ = new SubCans[num_cover];
  }

  ~FiPEIndex() {
    auto qnum = q_graph_->getVerticesCount();
    delete[] visited_u;
    delete[] visited_v;
    for (ui i = 0; i < qnum; i++) {
      for (ui j = 0; j < qnum; j++) {
        while (index_[i][j].size() > 1) {
          delete index_[i][j].back();
          index_[i][j].pop_back();
        }
      }
      delete[] index_[i];
    }
    delete[] index_;
    delete[] valid_cans_;
    delete embedding;
    delete indepInfo;
    delete[] subInfo_;
    delete[] subCans_;
  }

  // get valid cans of u_2 when select v for u_1, 24-11-7
  vector<VertexID> getNeighbors(VertexID u_1, VertexID u_2, VertexID v) {
#ifdef ANALYZE_TIME
    auto start_time = TimeOp::getClockNan();
#endif
    auto v_idx = b_search<vector<ui>>::search(valid_cans_[u_1].front(), valid_cans_[u_1].front().size(), v);
    if (v_idx == (ui)-1) {
      cout << "can't find " << v << " in valid_cans[" << u_1 << "]" << endl;
      exit(1);
    }
    auto &edges = *(index_[u_1][u_2].back());
    auto nbrs = edges.edge_ + edges.offset_[v_idx];
    auto nbrs_cnt = edges.offset_[v_idx + 1] - edges.offset_[v_idx];

    auto &old_cans = valid_cans_[u_2].back();

    // compute remained vertices
    auto result = std::move(SetOp::intersectTwo(old_cans, nbrs, nbrs_cnt));
#ifdef ANALYZE_TIME
    getNeighbors_time += TimeOp::getClockNan() - start_time;
#endif
    return result;
  }
};

#endif
