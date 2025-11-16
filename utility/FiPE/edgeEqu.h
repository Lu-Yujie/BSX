// edge equivalent information
#ifndef EDGE_EQUIVALENT_H
#define EDGE_EQUIVALENT_H
#include "graph/graph.h"
#include "bsx/SetOp.h"
#include<unordered_map>
using namespace std;

/**
 * If the current number of candidates is too large,
 * randomly select a subset of points from the current candidates
 */
class SubCans {
public:
  bool splitted;
  ui cur_start;
  ui cur_end;
  bool up_changed;
};

/**
 * store edge substitutable info of u and nxt
 * two sets of structures: cans_* filled at u, sub_* filled at nxt
 *     compute edge group by use these two sets alternatively
 * cur_s: point to the current group, group_id = cur_s
 *
 * _cnt: #groups
 * _offset: range of up_cans of each group
 * _up_idxs: candidates idxs list of each group(initial 0,1,...,n)
 * _down_offset: range of down_cans of each up
 * _edge_idxs: idxs of edges(c_down, initial 0,1,...,n)
 */
class EdgeSub {
public:
  ui c_cnt;
  vector<ui> edge_up_offset;
  vector<ui> edge_up_idxs;
  vector<ui> edge_down_offset;
  vector<ui> edge_down_idxs;
  //indicates edge_down is recorded to down_cans or not, store idx of v's
  static vector<ui> down_record;  // idx of v in data graph -> idx of v in down_cans

  ui p_cnt;
  vector<ui> p_edge_up_offset;
  vector<ui> p_edge_up_idxs;
  vector<ui> p_edge_down_offset;
  vector<ui> p_edge_down_idxs;

  vector<ui> down_group_num;  // group_id of each down_cans
  vector<ui> up_group_num;  // group_id of each up_cans

  ui cur_s;
};

/**
 * all kinds of unmatched nbrs
 * up_: all nbrs of up
 * down_: all nbrs of down
 * shared_: shared nbrs of u & nxt, compute the valid_cans of them based on up&down
 * up_indep_: up_indep_nbrs - shared_nbrs, only used for first edge, aka, up = start_vertex
 *      valid_cans based on up_cans
 * delayed_: if connected: nxt if (nxt not connected to up) else none
 *           else: (down_nbrs ∩ d_nbrs) - up_nbrs
 *           valid_cans: union of current group of down_cans
 *           if nxt is connected to up&down, valid_cans of nxt do not need union op.
 * down_indep_: down_nbrs - shared_nbrs - delayed_nbrs
 *              valid_cans: based on down_cans
 */
class Nbrs{
public:
  vector<VertexID> up_indep_;
  vector<VertexID> shared_;
  vector<VertexID> down_indep_;
  vector<VertexID> delayed_;

  void print() {
    cout << "up_indep_: ";
    for (auto&com:up_indep_) cout << com << ", ";
    cout << "\nshared_: ";
    for (auto&com:shared_) cout << com << ", ";
    cout << "\ndown_indep_: ";
    for (auto&com:down_indep_) cout << com << ", ";
    cout << "\ndelayed_: ";
    for (auto&com:delayed_) cout << com << ", ";
    cout << endl;
  }
};

/**
 * subs: substitutable info of each edge
 *       use can_idx as index, because influenced is based on the order of sub_cans
 * edges: record the matched edges: v_idx∈idx of cans of cur_s(u-1) -> v∈cans of cur_s(u)
 * when backtracking to u_n, re-compute (the sub_cans of u_n+1) & (c(u_n+1) on sub_cans)
 *
 * up2down: candidates list of each up to down
 * down_idxs: indicates idx of edge_down in down_cans
 * connected: up&down connected
 * same_nxt: if nxt in com_list, sub_cans of nxt do not need union operation
 * com_list: up_indep_(if start_vertex) + shared_ + down_indep_
 *
 * delayed_nbrs: delayed_nbrs = (nbrs_u - com_list) ∩ unmatched_u
 * delayed_inf: delayed_nbr is influenced or not in the current space
 *
 * the length of grouped, inf, inf_cans is equal to the number of edges
 * grouped: indicates whether the can is grouped
 * influenced: influenced u_nbrs of each edges, free at the end, [u_nbr][down_idx]
 * influenced_cans: cans of influenced_u
 */
class SubInfo {
public:
  EdgeSub subs;
  vector<VertexID> up_cans;
  vector<VertexID> down_cans;
  vector<VertexID> up2down;
  vector<ui> up_idxs;  // edge_idx -> up_can_idx
  vector<ui> down_idxs;  // edge_idx -> down_can_idx
  vector<bool> down_valid;
  vector<bool> up_valid;
  vector<bool> edge_valid;
  // auxiliary variables for compute sub info
  bool connected;
  Nbrs nbrs;

  vector<bool> grouped;
  vector<ui> nxt_up_idxs;  // 
  vector<ui> last_down_idxs;  // 
  vector<bool> delayed_inf;
  vector<bool>* influenced;
  vector<vector<VertexID>>* inf_cans;
  // valid_cans of shared nbrs of (up, down) respectively
  vector<vector<VertexID>>* up_shared_valid_cans;
  vector<vector<VertexID>>* down_shared_valid_cans;

  // compute com_list
  void init(const Graph* g, VertexID* order, ui depth, bool* visited_u, ui num_cover) {
    auto& up = order[depth];
    auto& down = order[depth+1];
    ui u_nbrs_cnt;
    auto u_nbrs = g->getVertexNeighbors(up, u_nbrs_cnt);
    ui d_nbrs_cnt;
    auto d_nbrs = g->getVertexNeighbors(down, d_nbrs_cnt);
    vector<VertexID> u_nbrs_unmatched;
    vector<VertexID> d_nbrs_unmatched;
    for (ui i = 0; i < u_nbrs_cnt; i++) {
      if (!visited_u[u_nbrs[i]]) u_nbrs_unmatched.emplace_back(u_nbrs[i]);
    }
    for (ui i = 0; i < d_nbrs_cnt; i++) {
      if (!visited_u[d_nbrs[i]]) d_nbrs_unmatched.emplace_back(d_nbrs[i]);
    }
    nbrs.shared_ = std::move(SetOp::intersectTwo(u_nbrs_unmatched, d_nbrs_unmatched));
    if (depth == 0) {
      nbrs.up_indep_ = std::move(SetOp::setDifference(u_nbrs, u_nbrs_cnt, nbrs.shared_, true));
      ui valid_up_indep_idx = 0;
      for (ui i = 0 ; i < nbrs.up_indep_.size(); i++)
        if (nbrs.up_indep_[i] != down)nbrs.up_indep_[valid_up_indep_idx++] = nbrs.up_indep_[i];
      nbrs.up_indep_.resize(valid_up_indep_idx);
    }
    if (g->checkEdgeExistence(up, down)) connected = true;
    else connected = false;

    if (depth + 2 >= num_cover) {  // means last edge
      nbrs.down_indep_.swap(d_nbrs_unmatched);
      SetOp::setDifference(nbrs.down_indep_, nbrs.shared_);
      // no nbrs.delayed_
    } else {
      auto nxt = order[depth+2];
      if (g->checkEdgeExistence(down, nxt)) {  // delay: nxt if up&nxt not connected
        if (!g->checkEdgeExistence(up,nxt)) {
          nbrs.delayed_.emplace_back(nxt);
        }
      } else {  // delay: (down_nbrs ∩ nxt_nbrs) - up_nbrs
        ui nxt_nbrs_cnt;
        auto nxt_nbrs = g->getVertexNeighbors(nxt, nxt_nbrs_cnt);
        nbrs.delayed_ = std::move(SetOp::intersectTwo(d_nbrs_unmatched, nxt_nbrs, nxt_nbrs_cnt));
        SetOp::setDifference(nbrs.delayed_, u_nbrs, u_nbrs_cnt);
      }
      // down_indep_: down_nbrs - shared_nbrs - delayed_nbrs
      nbrs.down_indep_.swap(d_nbrs_unmatched);
      SetOp::setDifference(nbrs.down_indep_, nbrs.shared_);
      SetOp::setDifference(nbrs.down_indep_, nbrs.delayed_);
    }
    // nbrs.print();

    // init influenced
    auto vnum = g->getVerticesCount();
    influenced = new vector<bool>[vnum];
    inf_cans = new vector<vector<VertexID>>[vnum];
    up_shared_valid_cans = new vector<vector<VertexID>>[vnum];
    down_shared_valid_cans = new vector<vector<VertexID>>[vnum];
  }

  ~SubInfo() {
    delete[] influenced;
    delete[] inf_cans;
    delete[] up_shared_valid_cans;
    delete[] down_shared_valid_cans;
  }
};

#endif
