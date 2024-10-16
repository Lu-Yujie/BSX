/**
 * compute node similarity
 * use LSH idea, compute lsh value based on adj vector
 * use sparse vector, the nbrs lists are just the sparse manifest. of adj vector
 * use landmarks, compute value for each landmark
 */
#ifndef NODE_SIMILARITY_H
#define NODE_SIMILARITY_H
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include "bsx.h"
using namespace std;
namespace NodeSim {
  #define MAX_LANDMARKS 10.0
  ///////////////////////// one node has one group nbrs(single vector) ////////////////////////////
  /**
   * choose landmarks based on nbrs_cnt(#log2(n), at most MAX_LANDMARKS)
   * sorted on nbrs_cnt, then equidistant sampling
   * store idx instead of vid
   */
  inline vector<ui>
  chooseLandmark(const ui *nbrs_cnt, ui num_node) {
    vector<ui> landmarks;
    double log10_2 = log10(2.0);
    double max_num = log10(num_node)/log10_2;
    double gap = num_node/(min(max_num, MAX_LANDMARKS));
    landmarks.reserve(gap);
    // sort the nbrs_cnt
    vector<pair<ui, ui>> nbrs_sorted;
    nbrs_sorted.reserve(num_node);
    for (ui i = 0; i < num_node; i++) {
      nbrs_sorted.emplace_back(nbrs_cnt[i], i);
    }
    sort(nbrs_sorted.begin(), nbrs_sorted.end(), [](pair<ui, ui> l, pair<ui, ui> r) {
      return l.first > r.first;
    });
    for (ui i = 0; i < num_node; i += gap) {
      ui landmark_idx = nbrs_sorted[i].second;
      landmarks.emplace_back(landmark_idx);
    }
    return move(landmarks);
  }
  /**
   * compute the similarity of each landmark(max:9):lms0,lsm1,...,lsm15(landmark score)
   * final score: lms0*10^15+lms1*10^14+...+lms14*10^1+lms15*10^0
   * note that the nbrs are in ascending order by default.
  */
  vector<int64_t>
  computeLSH(const ui *nbrs_cnt, const VertexID **nbrs, ui num_node, vector<ui> landmarks) {
    vector<int64_t> lsh_values;
    lsh_values.reserve(num_node);
    for (ui i = 0; i < num_node; i++) {
      int64_t score = 0;
      int scale = 1e15;
      for (auto& landmark_idx:landmarks) {
        // compute the #same_nbrs for each landmarks
        int lscore = 0;
        if (landmark_idx == i) lscore = nbrs_cnt[i];
        else {
          ui i_idx = 0, l_idx = 0;
          while (i_idx < nbrs_cnt[i] && l_idx < nbrs_cnt[landmark_idx] && lscore < 9) {
            if (nbrs[i][i_idx] == nbrs[landmark_idx][l_idx]) {
              lscore++;
              i_idx++;
              l_idx++;
            } else if (nbrs[i][i_idx] < nbrs[landmark_idx][l_idx]) {
              i_idx++;
            } else {
              l_idx++;
            }
          }
        }
        score += scale * lscore;
        scale /= 10;
      }
      lsh_values.emplace_back(score);
    }
    return move(lsh_values);
  }

  vector<int64_t>
  nodeSim(const ui *nbrs_cnt, const VertexID **nbrs, ui num_node) {
    // compute landmarks
    vector<ui> landmarks = move(chooseLandmark(nbrs_cnt, num_node));
    // compute the lsh value
    return move(computeLSH(nbrs_cnt, nbrs, num_node, landmarks));
  }

  //////////////////// single node has multi group nbrs(nested vector) ////////////////////////////
  //////////////////// adapted to BSXIndex structure                  ////////////////////////////

  vector<int64_t>  // compute node similarity of valid_cans of u
  nodeSim(BSXIndex& index, VertexID u) {
    // compute total #nbrs of each node
    auto& num_node = index.valid_cnt_[u].top();
    ui * nbrs_cnt = new ui[num_node];
    memset(nbrs_cnt, 0, sizeof(ui)*num_node);
    ui unbrs_cnt;
    auto unbrs = index.q_graph_->getVertexNeighbors(u, unbrs_cnt);
    for (ui v = 0; v < num_node; v++) {
      for (ui i = 0; i < unbrs_cnt; i++) {
        auto& unbr = unbrs[i];
        auto& edges = *(index.index_[u][unbr].top());
        nbrs_cnt[v] += edges.offset_[v+1] - edges.offset_[v];
      }
    }
    // compute landmarks
    vector<ui> landmarks = move(chooseLandmark(nbrs_cnt, num_node));
    cout << "landmarks: ";
    for (auto&tmp_ele:landmarks) cout << tmp_ele << ", ";
    cout << endl;
    // compute the lsh value
    vector<int64_t> lsh_values;
    lsh_values.reserve(num_node);
    for (ui i = 0; i < num_node; i++) {
      int64_t score = 0;
      int scale = exp10(MAX_LANDMARKS);
      for (auto& landmark_idx:landmarks) {
        // compute the #same_nbrs for each landmarks
        int lscore = 0;
        for (ui unbrs_idx = 0; unbrs_idx < unbrs_cnt; unbrs_idx++) {
          auto& unbr = unbrs[unbrs_idx];
          auto& edges = *(index.index_[u][unbr].top());
          auto i_nbrs = edges.edge_ + edges.offset_[i];
          auto i_cnt = edges.offset_[i+1] - edges.offset_[i];
          auto landmark_nbrs = edges.edge_ + edges.offset_[landmark_idx];
          auto landmark_cnt = edges.offset_[landmark_idx+1] - edges.offset_[landmark_idx];
          if (landmark_idx == i) lscore = nbrs_cnt[i];
          else {
            ui i_idx = 0, l_idx = 0;
            while (i_idx < i_cnt && l_idx < landmark_cnt && lscore < 9) {
              if (i_nbrs[i_idx] == landmark_nbrs[l_idx]) {
                lscore++;
                i_idx++;
                l_idx++;
              } else if (i_nbrs[i_idx] < landmark_nbrs[l_idx]) {
                i_idx++;
              } else {
                l_idx++;
              }
            }
          }
        }
        score += scale * lscore;
        scale /= 10;
      }
      lsh_values.emplace_back(score);
    }
    return move(lsh_values);
  }

} // namespace NodeSimilarity

#endif  // NODE_SIMILARITY_H
