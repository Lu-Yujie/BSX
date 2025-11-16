#ifndef SRC_VEQ_H
#define SRC_VEQ_H

#include "include/run.h"
#include "graph/graph.h"

namespace VEQ {

/// @brief init data graph & timeLimit for veq
/// @param d_graph data graph
/// @param time_limit The unit is milliseconds
void init_veq(std::string d_graph, int64_t time_limit);

/// @brief find q_graph on given dataGraph
/// @param q_graph query graph
/// @param embedding_cnt necessary variable, total num generated in limited_time
/// @param time_total necessary variable, total time, the unit is second
/// @return overtime or not
bool src_veq(std::string q_graph, mpz_t embedding_cnt, int64_t& time_total, uint64_t output_limit);

}  // namespace VEQ

#endif  //SRC_VEQ_H
