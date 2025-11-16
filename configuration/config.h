

#ifndef SUBGRAPHMATCHING_CONFIG_H
#define SUBGRAPHMATCHING_CONFIG_H


/**
 * Set the maximum size of a query graph. By default, we set the value as 64.
 */
#define MAXIMUM_QUERY_GRAPH_SIZE 64
#define HASH_TABLE_RATIO 1.2

/**
 * Setting the value as 1 is to (1) enable the neighbor label frequency filter (i.e., NLF filter); and (2) enable
 * to check the existence of an edge with the label information. The cost is to (1) build an unordered_map for each
 * vertex to store the frequency of the labels of its neighbor; and (2) build the label neighbor offset.
 * If the memory can hold the extra memory cost, then enable this feature to boost the performance. Otherwise, disable
 * it by setting this value as 0.
 */
#define OPTIMIZED_VLABELED_GRAPH 1

/**
 * have edge label or not
*/
// #define ELABELED_GRAPH

/**
 * Set intersection method.
 * 0: Hybrid method; 1: Merge based set intersections.
 */
#define HYBRID 0

/**
 * Accelerate set intersection with SIMD instructions.
 * 0: AVX2; 1: AVX512; 2: Basic;
 */
#define SI 0

/**
 * Define ENABLE_FAILING_SET to enable the failing set pruning set intersection method.
 */
#define ENABLE_FAILING_SET

/**
 * Define DYNAMIC_FAIL_REASON to enable the failing set pruning set intersection method.
 */
// #define DYNAMIC_FAIL_REASON

/**
 * Define ENABLE_EQUIVALENT_SET to enable the equivalent set pruning technique in VEQ
 */
#define ENABLE_EQUIVALENT_SET

/**
 * Define BSX_SIM_THRESHOLD to change the similiarity threshold, 100 means identical.
 */
#define BSX_SIM_THRESHOLD 25

/**
 * Define PART_SEP_KIND to choose diff seperation methods.
 */
#define PART_SEP_KIND 1

#define PRINT_SEPARATOR "------------------------------"

#endif //SUBGRAPHMATCHING_CONFIG_H
