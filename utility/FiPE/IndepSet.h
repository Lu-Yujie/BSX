// compute min vertex cover based on the #cans
#ifndef FIPE_INDEP_H
#define FIPE_INDEP_H

#include "graph/graph.h"
#include "gmp.h"
using namespace std;

class FiPEIndep {
  ui* indep_con_cnt;    // count the number of times each indep_cans may conflict
  ui** sep_flag;        // seperate indep cans
  mpz_t label_embeddings;  // #embeddings of one kind of label, for enumeration
  bool* visited_v;

  // for memory free
  ui qnum_;
public:
  // Quickly determine indep
  bool* indep_bool;
  ui num_cover_;
  bool* used_cans_;
  vector<vector<VertexID>> cans;

  // get embedding result
  mpz_t embedding_step;   // #embeddings of one step
  mpz_t embedding_total;  // #embeddings of one enumeration
  mpz_t embedding_uncon;  // #embeddings of the step without potential conflict
  mpz_t remained;
  bool uncon;             // indicates embedding_uncon is search or not

  FiPEIndep(ui qnum, ui dnum, ui num_cover, VertexID* order) {
    qnum_ = qnum;
    num_cover_ = num_cover;
    cans.resize(qnum);
    visited_v = new bool[dnum];
    memset(visited_v, false, sizeof(bool)*dnum);
    indep_con_cnt = new ui[dnum];
    memset(indep_con_cnt, 0, sizeof(ui)*dnum);
    indep_bool = new bool[qnum];
    memset(indep_bool, true, sizeof(bool)*qnum);
    for (ui i = 0; i < num_cover; i++) {
      cans[order[i]].emplace_back(0);
      indep_bool[order[i]] = false;
    }
    used_cans_ = new bool[dnum];
    sep_flag = new ui*[qnum];
    for (ui i = 0; i < qnum; i++) {
      // each node has 3 seperation point, seperate into four parts
      // 1 upward conflict: seperate nodes based on whether it shows upward
      // 2 downward   ''  :    ''     ''     ''        ''       ''   downward
      // 4 parts: up-down,up-x,x-down,x-x
      sep_flag[i] = new ui[3];
    }
    mpz_init(embedding_step);
    mpz_init(embedding_total);
    mpz_init(embedding_uncon);
    mpz_init(label_embeddings);
    mpz_init(remained);
  }
  ~FiPEIndep() {
    delete[] indep_con_cnt;
    delete[] indep_bool;
    delete[] used_cans_;
    for (ui i = 0; i < qnum_; i++) {
      delete[] sep_flag[i];
    }
    delete[] sep_flag;
    delete[] visited_v;
    mpz_clear(embedding_step);
    mpz_clear(embedding_total);
    mpz_clear(embedding_uncon);
    mpz_clear(label_embeddings);
    mpz_clear(remained);
  }

  void homoEnum() {
    mpz_set_ui(embedding_step, 1);
    for (auto can:cans) {
      mpz_mul_ui(embedding_step, embedding_step, can.size());
    }
  }

  void enumeration(const Graph*q_graph) {
    auto label_num = q_graph->getLabelsCount();
    mpz_set_ui(embedding_step, 1);
    for (ui l_idx = 0; l_idx < label_num; l_idx++) {
        ui nodes_num;
        // these nodes have the same label
        auto nodes = q_graph->getVerticesByLabel(l_idx, nodes_num);
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
            mpz_mul_ui(embedding_step, embedding_step, cans[nodes[0]].size());
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
            for (ui j = 0; j < cans[node].size(); j++) indep_con_cnt[v_cans[j]]++;
        }
        // 2.second scan, compute downward conflict
        for (ui i = 0; i < nodes_num; i++) {
            auto& node = nodes[i];
            auto& v_cans = cans[node];
            auto v_cans_cnt = cans[node].size();
            auto& downward_sep0 = sep_flag[node][0];  // 0->sep the up-conflicts
            // auto& downward_sep1 = sep_flag[cur_idx][2];  // 1->sep the up-uncon.
            int forward_idx = 0;
            // int middle_idx = sep_flag[cur_idx][1];
            // int backward_idx = v_cans_cnt - 1;
            for (ui j = 0; j < v_cans_cnt; j++) indep_con_cnt[v_cans[j]]--;
            downward_sep0 = sepDiff(v_cans, indep_con_cnt, forward_idx, v_cans_cnt - 1);
        }
        // 3.enumerate the nodes based on diff features of 4 parts
        // ** just 2 parts so far
        enum4Parts(sep_flag, nodes, nodes_num, cans, visited_v, label_embeddings);
        mpz_mul(embedding_step, embedding_step, label_embeddings);
        if (mpz_cmp_ui(embedding_step, 0) == 0) return;
    }
    return;
  }

  // according to indep_con_cnt info, seperate v_cans into two parts, return the #first_part(true)
  static ui sepDiff(std::vector<VertexID> &v_cans, const ui *indep_con_cnt,
                    int forward_idx, int backward_idx) {
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

  /**
   * 4 parts: up-down,up-x,x-down,x-x; down&x 2 parts so far
   * TODO: opt to three parts
   */
  static void enum4Parts(ui **&sep_flags, const VertexID* nodes, ui nodes_num,
                        std::vector<std::vector<VertexID>>& cans, bool *&visited_v, mpz_t cur_cnt) {
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

  // select indepent nodes according to #cans
  static ui indepSetOnCans(const Graph* graph, VertexID *order, ui* cans_cnt) {
    ui cover_num = 0;
    ui n = graph->getVerticesCount();
    auto offset = graph->getOffsets();
    auto edges = graph->getEdges();
    ui *degree = new ui[n];
    bool *visited = new bool[n];
    memset(visited, false, sizeof(bool)*n);

    for (VertexID i = 0; i < n; i++) {
      degree[i] = graph->getVertexDegree(i);
    }

    // select (unvisited && degree > 0) node with least #cans, then max degree
    while (true) {
      ui maxDegree = 0, minCans = (ui)-1;
      VertexID selectedNode = (VertexID)-1;
      for (VertexID i = 0; i < n; i++) {
        if (!visited[i] && degree[i] > 0) {
          if (cans_cnt[i] < minCans
              || (cans_cnt[i] == minCans && maxDegree < degree[i])) {
            minCans = cans_cnt[i];
            maxDegree = degree[i];
            selectedNode = i;
          }
        }
      }

      if (selectedNode == (VertexID)-1) {
        break;  // all edges are covered
      }

      visited[selectedNode] = true;
      order[cover_num++] = selectedNode;

      // delete joint edges
      for (ui i = offset[selectedNode]; i < offset[selectedNode + 1]; i++) {
        VertexID neighbor = edges[i];
        degree[neighbor]--;
      }
    }

    if (cover_num == 1) {
      ui nbrs_cnt;
      auto nbrs = graph->getVertexNeighbors(order[0], nbrs_cnt);
      order[cover_num++] = nbrs[0];
    }

    delete[] degree;
    delete[] visited;
    return cover_num;
  }

  // select indepent nodes according to edge connection, for edge substitutable
  // 1st_u: min degree, min cans
  // nth_u: connect to (n-1)th, min degree, min cans
  static ui indepSetOnEdge(const Graph* graph, VertexID *order, ui* cans_cnt) {
    ui cover_num = 0;
    ui n = graph->getVerticesCount();
    auto offset = graph->getOffsets();
    auto edges = graph->getEdges();
    vector<int> degree(n, 0);
    bool *visited = new bool[n];
    memset(visited, false, sizeof(bool)*n);
    bool *connected = new bool[n];
    memset(connected, false, sizeof(bool)*n);

    for (VertexID i = 0; i < n; i++) {
      if (graph->getVertexDegree(i) > 1)
        degree[i] = graph->getVertexDegree(i);
    }

    // select (unvisited && degree > 0) node with least #cans, then max degree
    while (true) {
      ui minDegree = (ui)-1, minCans = (ui)-1, coned = false;
      VertexID selectedNode = (VertexID)-1;
      for (VertexID i = 0; i < n; i++) {
        if (!visited[i] && degree[i] > 0) {
          if ((connected[i] && !coned)  // 1.connected better
              || ((connected[i] == coned)  // 2.same connection
                  && (minDegree > degree[i]  // 2.1 small degree better
                      || (minDegree == degree[i]  // 2.2 same degree
                          && cans_cnt[i] < minCans)))) {  // 2.2.1 small #cans better
            minCans = cans_cnt[i];
            minDegree = degree[i];
            selectedNode = i;
            coned = connected[i];
          }
        }
      }

      if (selectedNode == (VertexID)-1) {
        break;  // all edges are covered
      }

      visited[selectedNode] = true;
      order[cover_num++] = selectedNode;
      memset(connected, false, sizeof(bool)*n);

      // delete joint edges
      for (ui i = offset[selectedNode]; i < offset[selectedNode + 1]; i++) {
        VertexID neighbor = edges[i];
        degree[neighbor]--;
        connected[neighbor] = true;
      }
    }

    if (cover_num == 1) {
      ui nbrs_cnt;
      auto nbrs = graph->getVertexNeighbors(order[0], nbrs_cnt);
      order[cover_num++] = nbrs[0];
    }

    delete[] visited;
    delete[] connected;
    return cover_num;
  }

  // select indepent nodes according to degree
  static ui indepSetOnDegree(const Graph* graph, VertexID *order) {
    ui cover_num = 0;
    ui n = graph->getVerticesCount();
    auto offset = graph->getOffsets();
    auto edges = graph->getEdges();
    vector<int> degree(n, 0);
    vector<bool> visited(n, false);

    for (VertexID i = 0; i < n; i++) {
      if (graph->getVertexDegree(i) > 1)
        degree[i] = graph->getVertexDegree(i);
    }

    // select node with largest degree
    while (true) {
      int maxDegree = 0, selectedNode = (VertexID)-1;
      for (VertexID i = 0; i < n; i++) {
        if (!visited[i] && degree[i] > maxDegree) {
          maxDegree = degree[i];
          selectedNode = i;
        }
      }

      if (selectedNode == (VertexID)-1) {
          break;  // all edges are covered
      }

      visited[selectedNode] = true;
      order[cover_num++] = selectedNode;

      // delete joint edges
      for (ui i = offset[selectedNode]; i < offset[selectedNode + 1]; i++) {
        VertexID neighbor = edges[i];
        degree[neighbor]--;
      }
    }

    if (cover_num == 1) {
      ui nbrs_cnt;
      auto nbrs = graph->getVertexNeighbors(order[0], nbrs_cnt);
      order[cover_num++] = nbrs[0];
    }

    return cover_num;
  }

};

#endif