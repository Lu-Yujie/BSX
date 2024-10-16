/**
 * compute the max indepedent set
 * code from https://github.com/LijunChang/Near-Maximum-Independent-Set.
 * Lijun Chang, Wei Li and Wenjie Zhang.
 * Computing A Near-Maximum Independent Set in Linear Time by Reducing-Peeling.
 * changes:
 *  1. support graph format
 *  2. return specific set contents
 *  3. remove unused code
*/

#ifndef IndepSet_H
#define IndepSet_H
#include "graph/graph.h"

using namespace std;

class IndepSet {
  ui n, m; //number of nodes and edges of the graph
	const ui *pstart; //offset of neighbors of nodes
	VertexID *edges; //adjacent ids of edges

public:
  IndepSet(const Graph* graph) {
    n = graph->getVerticesCount();
    m = graph->getEdgesCount() * 2;
    pstart = graph->getOffsets();
    edges = new VertexID[m];
    copy(graph->getEdges(), graph->getEdges()+m, edges);
  }
  ~IndepSet() {
    delete[] edges;
  }

private:
	bool check_is(const bool *is, int count) {
    cout << "result: ";
    for(ui i = 0; i < n; i++) cout << is[i] << ", ";
    cout << endl;
    int cnt = 0;
    for(ui i = 0;i < n;i ++) if(is[i]) ++ cnt;
    if(count != cnt) {
      printf("WA count error! %d\n", cnt);
      return false;
    }

    bool maximal = true;
    for(ui i = 0;i < n;i ++) {
      if(is[i]) {
        for(ui j = pstart[i];j < pstart[i+1];j ++) if(is[edges[j]]) {
          cout << "WA conflict on " << i << ", " << edges[j] << endl;
          return false;
        }
      }
      else if(maximal) {
        bool find = false;
        for(ui j = pstart[i];j < pstart[i+1];j ++) if(is[edges[j]]) {
          find = true;
          break;
        }
        if(!find) {
          maximal = false;
        }
      }
    }
    if(!maximal) {
      printf("** WA is not maximal\n");
      return false;
    }
    else printf("WA is a maximal\n");
    return true;
  }

  int delete_vertex(ui v, const ui *pend, bool *is, int *degree, vector<ui> &degree_ones, vector<ui> &degree_twos) {
    is[v] = false;
    int res = 0;
    for(int k = pstart[v];k < pend[v];k ++) if(is[edges[k]]) {
      int w = edges[k];
      -- degree[w];
      if(degree[w] == 0) ++ res;
      else if(degree[w] == 1) degree_ones.push_back(w);
      else if(degree[w] == 2) degree_twos.push_back(w);
    }
    return res;
  }

  int exist_edge(ui u, ui v, const ui *pend) {
    if(pend[u]-pstart[u] < pend[v]-pstart[v]) {
      for(ui i = pstart[u];i < pend[u];i ++) {
        if(edges[i] == v) return 1;
      }
      return 0;
    }
    for(ui i = pstart[v];i < pend[v];i ++) {
      if(edges[i] == u) return 1;
    }
    return 0;
  }

  ui edge_rewire(ui u, const ui *pend, ui v, ui w) {
  #ifndef NDEBUG
    for(ui i = pstart[u];i < pend[u];i ++) if(edges[i] == w) printf("WA preexist edge!\n");
  #endif

    for(ui i = pstart[u];i < pend[u];i ++) if(edges[i] == v) {
      edges[i] = w;
      return i;
    }
    printf("WA in edge_rewire!\n");
    return 0;
  }

  void shrink(ui u, ui &end, const bool *is) {
    ui i = pstart[u];
    while(true) {
      while(i < end&&is[edges[i]]) ++ i;
      while(i < end&&!is[edges[end-1]]) -- end;

      if(i >= end) break;
      swap(edges[i], edges[end-1]);
    }
  }

public:
  // clear returned bool*, true->indep, false->o.w.
  pair<bool*, ui> linearTime() {
    bool *is = new bool[n];
    for(ui i = 0;i < n;i ++) is[i] = true;

    int *bin_head = new int[n];
    int *bin_next = new int[n];
    int *degree = new int[n];
    memset(bin_head, -1, sizeof(int)*n);

    vector<ui> degree_ones, degree_twos;
    vector<pair<ui,ui>> S;
    vector<pair<pair<ui,ui>, ui>> modified_edges;

    int max_d = 0, res = 0;
    for(ui i = 0;i < n;i ++) {
      degree[i] = pstart[i+1] - pstart[i];
      bin_next[i] = bin_head[degree[i]];
      bin_head[degree[i]] = i;

      if(degree[i] == 0) ++ res;
      else if(degree[i] == 1) degree_ones.push_back(i);
      else if(degree[i] == 2) degree_twos.push_back(i);

      if(degree[i] > max_d) max_d = degree[i];
    }

    char *fixed = new char[n];
    memset(fixed, 0, sizeof(char)*n);

    ui *pend = new ui[n];
    for(ui i = 0;i < n;i ++) pend[i] = pstart[i+1];

    int kernal_size = 0, inexact = 0, first_time = 1, S_size = (int)S.size();
    int kernal_edges = 0;
    while(!degree_ones.empty()||!degree_twos.empty()||max_d >= 3) {
      while(!degree_ones.empty()||!degree_twos.empty()) {
        while(!degree_ones.empty()) {
          ui u = degree_ones.back();
          degree_ones.pop_back();
          if(!is[u]||degree[u] != 1) continue;

          int cnt = 0;
          for(int j = pstart[u];j < pend[u];j ++) if(is[edges[j]]) {
            ++ cnt;
            res += delete_vertex(edges[j], pend, is, degree, degree_ones, degree_twos);
          }
        }

        while(!degree_twos.empty()&&degree_ones.empty()) {
          ui u = degree_twos.back();
          degree_twos.pop_back();
          if(!is[u]||degree[u] != 2) continue;

          shrink(u, pend[u], is);
          ui u1 = edges[pstart[u]], u2 = edges[pstart[u]+1];

          ui pre = u, cnt = 1;
          while(u1 != u&&degree[u1] == 2) {
            ++ cnt;
            shrink(u1, pend[u1], is);
            int tmp = u1;
            if(edges[pstart[u1]] != pre) u1 = edges[pstart[u1]];
            else u1 = edges[pstart[u1]+1];
            pre = tmp;
          }
          if(u1 == u) {
            res += delete_vertex(u, pend, is, degree, degree_ones, degree_twos);
            continue;
          }

          pre = u;
          while(degree[u2] == 2) {
            ++ cnt;
            shrink(u2, pend[u2], is);
            int tmp = u2;
            if(edges[pstart[u2]] != pre) u2 = edges[pstart[u2]];
            else u2 = edges[pstart[u2]+1];
            pre = tmp;
          }
          if(u1 == u2) {
            res += delete_vertex(u1, pend, is, degree, degree_ones, degree_twos);
            continue;
          }

          shrink(u1, pend[u1], is);
          shrink(u2, pend[u2], is);

          if(cnt%2 == 1) {
            if(exist_edge(u1, u2, pend)) {
              res += delete_vertex(u1, pend, is, degree, degree_ones, degree_twos);
              res += delete_vertex(u2, pend, is, degree, degree_ones, degree_twos);
            }
            else if(cnt > 1) {
              ui idx = pstart[pre];
              if(edges[idx] == u2) ++ idx;
              u = edges[idx];
              edges[idx] = u1;
              if(!first_time) modified_edges.push_back(make_pair(make_pair(pre,u), u1));

              u2 = pre;
              while(u != u1) {
                is[u] = 0;
                ui tmp = u;
                if(edges[pstart[u]] == pre) u = edges[pstart[u]+1];
                else u = edges[pstart[u]];
                S.push_back(make_pair(tmp, u));
                pre = tmp;
              }

              edge_rewire(u1, pend, pre, u2);
              if(!first_time) modified_edges.push_back(make_pair(make_pair(u1,pre), u2));
            }
          }
          else {
            ui v2 = pre, v1 = pre;
            pre = u2;
            while(v1 != u1) {
              is[v1] = false;
              int tmp = v1;
              if(edges[pstart[v1]] == pre) v1 = edges[pstart[v1]+1];
              else v1 = edges[pstart[v1]];
              S.push_back(make_pair(tmp, v1));
              pre = tmp;
            }
            v1 = pre;
            if(exist_edge(u1, u2, pend)) {
              -- degree[u1];
              -- degree[u2];
              if(degree[u1] == 2) degree_twos.push_back(u1);
              if(degree[u2] == 2) degree_twos.push_back(u2);
            }
            else {
              edge_rewire(u1, pend, v1, u2);
              edge_rewire(u2, pend, v2, u1);
              if(!first_time) {
                modified_edges.push_back(make_pair(make_pair(u1,v1), u2));
                modified_edges.push_back(make_pair(make_pair(u2,v2), u1));
              }
            }
          }
        }
      }

      if(first_time) {
        S_size = (int)S.size();
        first_time = 0;
        for(ui k = 0;k < n;k ++) {
          if(is[k]&&degree[k] > 0) {
            ++ kernal_size;
            for(ui j = pstart[k];j < pend[k];j ++) if(is[edges[k]]) ++ kernal_edges;
          }
          else fixed[k] = 1;
        }
      }

      while(degree_ones.empty()&&degree_twos.empty()) {
        while(max_d >= 3&&bin_head[max_d] == -1) -- max_d;
        if(max_d < 3) break;

        int v = -1;
        for(v = bin_head[max_d];v != -1;) {
          int tmp = bin_next[v];
          if(is[v]&&degree[v] > 0) {
            if(degree[v] < max_d) {
              bin_next[v] = bin_head[degree[v]];
              bin_head[degree[v]] = v;
            }
            else {
              S.push_back(make_pair(v,n)); ++ inexact;

              res += delete_vertex(v, pend, is, degree, degree_ones, degree_twos);

              bin_head[max_d] = tmp;
              break;
            }
          }
          v = tmp;
        }
        if(v == -1) bin_head[max_d] = -1;
      }
    }

    for(int i = S.size()-1;i >= 0;i --) {
      ui u1 = S[i].first, u2 = S[i].second;

      if(u2 != n) {
        if(!is[u2]) {
          is[u1] = true;
          ++ res;
        }
        continue;
      }

      int ok = 1;
      for(ui i = pstart[u1];i < pstart[u1+1];i ++) if(is[edges[i]]) {
        ok = 0;
        break;
      }
      if(ok) {
        is[u1] = true;
        ++ res;
      }
    }

    delete[] bin_head;
    delete[] bin_next;
    delete[] degree;
    delete[] pend;

    for(int i = (int)modified_edges.size()-1;i >= 0;i --) {
      ui u = modified_edges[i].first.first, u1 = modified_edges[i].first.second, u2 = modified_edges[i].second;
      for(ui j = pstart[u];j < pstart[u+1];j ++) if(edges[j] == u2) {
        edges[j] = u1;
        break;
      }
    }
    delete[] fixed;
    return make_pair(is, res);
  }

};

#endif
