#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "matchingcommand.h"
#include "graph/graph.h"
#include "GenerateFilteringPlan.h"
#include "FilterVertices.h"
#include "BuildEdgeIndex.h"
#include "GenerateQueryPlan.h"
#include "EvaluateQuery.h"
#include "sourceCode/VEQ/veq.h"
#include "pretty_print.h"

using namespace std;

/**
 * declare the sets to control the methods
 * 1. supported is the all available 'methods'
 * 2. incompatible set records all the incompatible 'combinations'
 * 3. wanted set records the 'combination' we want, superior to exclude set, inf. to incom.
 * 4. excluded record the 'methods' we don't care
 * 5. engine_wo_order: engine methods with dynamic order
 * 6. engine_wo_filter: engine methods without filter
 * 7. src_method: methods with source code
*/
class Methods {
public:
  enum Midx {
    E_IDX=0,
    F_IDX=1,
    O_IDX=2,
    S_IDX=3,
  };
  // source methods idx
  // vector for scanning, set for searching
  vector<vector<string>> supported;  // 0->engine, 1->filter, 2->order
  set<string> incompatible_set;
  set<string> wanted_set;
  vector<set<string>> excluded;
  set<string> engine_wo_order;
  set<string> engine_wo_filter;
  vector<string> src_method;
  vector<ui> cur_idx;  // used for scan the methods
  bool src_com;        // false->combanition, true->sourceCode

  /// @brief remove spaces at start and end of a string
  std::string trim(const std::string& str) {
    auto start = std::find_if_not(str.begin(), str.end(), [](unsigned char c) {
      return std::isspace(c);
    });
    auto end = std::find_if_not(str.rbegin(), str.rend(), [](unsigned char c) {
      return std::isspace(c);
    }).base();
    return (start < end ? std::string(start, end) : "");
  }

  /// @brief Function to split a string based on a delimiter and return a vector of tokens
  vector<string> split(const string& s, char delimiter) {
    vector<string> tokens;
    stringstream ss(s);
    string token;
    while (getline(ss, token, delimiter)) {
      tokens.push_back(trim(token));
    }
    return tokens;
  }

  string getCurMethods(Midx idx) {
    string engine;
    switch (idx) {
      case E_IDX:
        return supported[E_IDX][cur_idx[E_IDX]];
        break;
      case F_IDX:
        engine = supported[E_IDX][cur_idx[E_IDX]];
        return engine_wo_filter.find(engine) == engine_wo_filter.end()
               ? supported[F_IDX][cur_idx[F_IDX]] : "null";
        break;
      case O_IDX:
        engine = supported[E_IDX][cur_idx[E_IDX]];
        return engine_wo_order.find(engine) == engine_wo_order.end()
               ? supported[O_IDX][cur_idx[O_IDX]] : "null";
        break;
      case S_IDX:
        return src_method[cur_idx[S_IDX]];
        break;
      default:
        cerr << "wrong idx type" << endl;
        exit(-1);
        break;
    }
    return "null";
  }

  Methods (string conf_path) {
    src_com = false;
    auto supp_file = conf_path + "/support.txt";
    auto excluded_file = conf_path + "/excluded.txt";
    auto incompatible_file = conf_path + "/incompatible.txt";
    auto wanted_file = conf_path + "/wanted.txt";
    auto src_file = conf_path + "/src.txt";
    auto special_file = conf_path + "/special_engine.txt";
    supported.resize(3);
    excluded.resize(3);
    cur_idx.resize(4, 0);

    // read files and init variables
    string line;
    std::fstream methods_file;
    // 1. supported file
    methods_file.open(supp_file, std::ios::in);
    if (methods_file.is_open()) {
      while (getline(methods_file, line)) {
        if (line.empty())
          continue;
        uint64_t colonPos = line.find(':');
        if (colonPos != string::npos && colonPos + 1 < line.length()) {
          string contentAfterColon = line.substr(colonPos + 1); // Get content after colon
          vector<string> tokens = std::move(split(contentAfterColon, ',')); // Split by comma
          if (line.find("engine") != string::npos) {
            supported[E_IDX].insert(supported[E_IDX].end(), tokens.begin(), tokens.end());
          } else if (line.find("filter") != string::npos) {
            supported[F_IDX].insert(supported[F_IDX].end(), tokens.begin(), tokens.end());
          } else if (line.find("order") != string::npos) {
            supported[O_IDX].insert(supported[O_IDX].end(), tokens.begin(), tokens.end());
          }
        }
      }
    } else {
        cerr << "Error: Unable to open file " << supp_file << endl;
    }
    methods_file.close();
    // 2. excluded file
    methods_file.open(excluded_file, std::ios::in);
    if (methods_file.is_open()) {
      while (getline(methods_file, line)) {
        if (line.empty())
          continue;
        uint64_t colonPos = line.find(':');
        if (colonPos != string::npos && colonPos + 1 < line.length()) {
          string contentAfterColon = line.substr(colonPos + 1); // Get content after colon
          vector<string> tokens = std::move(split(contentAfterColon, ','));
          if (line.find("engine") != string::npos) {
            excluded[E_IDX].insert(tokens.begin(), tokens.end());
          } else if (line.find("filter") != string::npos) {
            excluded[F_IDX].insert(tokens.begin(), tokens.end());
          } else if (line.find("order") != string::npos) {
            excluded[O_IDX].insert(tokens.begin(), tokens.end());
          }
        }
      }
    } else {
        cerr << "Error: Unable to open file " << excluded_file << endl;
    }
    methods_file.close();
    // 3. incompatible file
    methods_file.open(incompatible_file, std::ios::in);
    if (methods_file.is_open()) {
      while (getline(methods_file, line)) {
        incompatible_set.insert(line);
      }
    } else {
        cerr << "Error: Unable to open file " << incompatible_file << endl;
    }
    methods_file.close();
    // 4. wanted file
    methods_file.open(wanted_file, std::ios::in);
    if (methods_file.is_open()) {
      while (getline(methods_file, line)) {
        wanted_set.insert(line);
      }
    } else {
        cerr << "Error: Unable to open file " << wanted_file << endl;
    }
    methods_file.close();
    // 5. src file
    methods_file.open(src_file, std::ios::in);
    if (methods_file.is_open()) {
      while (getline(methods_file, line)) {
        src_method.emplace_back(line);
      }
    } else {
        cerr << "Error: Unable to open file " << src_file << endl;
    }
    methods_file.close();
    // 6. special file
    methods_file.open(special_file, std::ios::in);
    if (methods_file.is_open()) {
      while (getline(methods_file, line)) {
        if (line.empty())
          continue;
        uint64_t colonPos = line.find(':');
        if (colonPos != string::npos && colonPos + 1 < line.length()) {
          string contentAfterColon = line.substr(colonPos + 1); // Get content after colon
          vector<string> tokens = std::move(split(contentAfterColon, ','));
          if (line.find("wo_order") != string::npos) {
            engine_wo_order.insert(tokens.begin(), tokens.end());
          } else if (line.find("wo_filter") != string::npos) {
            engine_wo_filter.insert(tokens.begin(), tokens.end());
          }
        }
      }
    } else {
        cerr << "Error: Unable to open file " << special_file << endl;
    }
    next(true);  // find the first method
  }

  /// @brief prepare for the next methods, add your special processing constraints here
  bool next(bool first = false) {
    string filter, order, engine;
    if (src_com) {
      cur_idx[S_IDX]++;
      return cur_idx[S_IDX] < src_method.size();
    }
    if (!first && cur_idx[E_IDX] < supported[E_IDX].size()) {  // if not first, find from the next available methods
      engine = supported[E_IDX][cur_idx[E_IDX]];
      // next combanition, based on cur_engine
      if (engine_wo_order.find(engine) != engine_wo_order.end()) {
        if (engine_wo_filter.find(engine) != engine_wo_filter.end()) {
          cur_idx[E_IDX]++;
        } else {
          if (cur_idx[F_IDX] == supported[F_IDX].size() - 1) {
            cur_idx[F_IDX] = 0;
            cur_idx[E_IDX]++;
          } else {
            cur_idx[F_IDX]++;
          }
        }
      } else {
        if (cur_idx[O_IDX] == supported[O_IDX].size() - 1) {
          if (engine_wo_filter.find(engine) != engine_wo_filter.end()) {
            cur_idx[E_IDX]++;
          } else {
            if (cur_idx[F_IDX] == supported[F_IDX].size() - 1) {
              cur_idx[F_IDX] = 0;
              cur_idx[E_IDX]++;
            } else {
              cur_idx[F_IDX]++;
            }
          }
          cur_idx[O_IDX] = 0;
        } else {
          cur_idx[O_IDX]++;
        }
      }
    }

    while (cur_idx[E_IDX] < supported[E_IDX].size()) {
      filter = supported[F_IDX][cur_idx[F_IDX]];
      order = supported[O_IDX][cur_idx[O_IDX]];
      engine = supported[E_IDX][cur_idx[E_IDX]];
      string combanition = filter + order + engine;
      bool flag_run = true;

      // check constrains
      do {
        if (incompatible_set.find(combanition) != incompatible_set.end() // incompatible methods
            || (engine == "RM" && engine != order)) {  // RM engine workes well only with order RM
            flag_run = false;
            break;
        }
        if (wanted_set.find(combanition) != wanted_set.end()) break;  // want this combination
        if ((engine_wo_filter.find(engine) == engine_wo_filter.end()
             && excluded[F_IDX].find(filter) != excluded[F_IDX].end())  // need filter & excluded
            || (engine_wo_order.find(engine) == engine_wo_order.end()
                && excluded[O_IDX].find(order) != excluded[O_IDX].end())  // need order & excluded
            || excluded[E_IDX].find(engine) != excluded[E_IDX].end()) {
            flag_run = false;
            break;
        }
      } while(0);

      if (flag_run) return true;

      // next combanition
      if (engine_wo_order.find(engine) != engine_wo_order.end()) {
        if (engine_wo_filter.find(engine) != engine_wo_filter.end()) {
          cur_idx[E_IDX]++;
        } else {
          if (cur_idx[F_IDX] == supported[F_IDX].size() - 1) {
            cur_idx[F_IDX] = 0;
            cur_idx[E_IDX]++;
          } else {
            cur_idx[F_IDX]++;
          }
        }
      } else {
        if (cur_idx[O_IDX] == supported[O_IDX].size() - 1) {
          if (engine_wo_filter.find(engine) != engine_wo_filter.end()) {
            cur_idx[E_IDX]++;
          } else {
            if (cur_idx[F_IDX] == supported[F_IDX].size() - 1) {
              cur_idx[F_IDX] = 0;
              cur_idx[E_IDX]++;
            } else {
              cur_idx[F_IDX]++;
            }
          }
          cur_idx[O_IDX] = 0;
        } else {
          cur_idx[O_IDX]++;
        }
      }
    }
    // no more combanition which is valid, start src
    src_com = true;
    cur_idx[S_IDX] = 0;
    return cur_idx[S_IDX] < src_method.size();
  }

  void print() {
    cout << "supported: " << supported << endl;
    cout << "excluded: " << excluded << endl;
    cout << "incompatible: " << incompatible_set.size() << endl;
    cout << "wanted: " << wanted_set << endl;
    cout << "src: " << src_method << endl;
    cout << "wo_filter: " << engine_wo_filter << endl;
    cout << "wo_order: " << engine_wo_order << endl;
  }
};

#endif  // EXPERIMENT_H
