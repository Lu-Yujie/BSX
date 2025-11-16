#include "veq.h"

namespace VEQ {

void init_veq(std::string d_graph, int64_t time_limit) {
  timeLimit = time_limit;
  mpq_init(nMatch);
  mpq_init(nCurrMatch);
  mpq_init(nRemainingMatch);
  mpq_init(nMaxMatch);
  mpq_set_d(nMaxMatch, numeric_limits<double>::max());
  mpq_init(mpqOne);
  mpq_set_d(mpqOne, 1);
  mpq_init(mpqTmp);
  mpq_init(mpqZero);
  mpq_set_d(mpqZero, 0);
  ReadIgraphFormat(d_graph, dataGraph);
  ProcessDataGraphs();
  AllocateForDataGraph();
}

bool src_veq(std::string q_graph, mpz_t embedding_cnt, int64_t& time_total, uint64_t output_limit) {
  // one query graph per time
  if (queryGraph.size() != 0) delete queryGraph.back();
  if (output_limit != numeric_limits<uint64_t>::max()) {
    mpq_set_d(nMaxMatch, output_limit);
  }
  queryGraph.clear();
  ReadIgraphFormat(q_graph, queryGraph);
  return ProcessQuery(embedding_cnt, time_total);
}

}  // namespace VEQ
