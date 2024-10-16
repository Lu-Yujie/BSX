#include <chrono>
#include <fstream>

#include "matchingcommand.h"
#include "graph/graph.h"
#include "GenerateFilteringPlan.h"
#include "FilterVertices.h"
#include "BuildEdgeIndex.h"
#include "GenerateQueryPlan.h"
#include "EvaluateQuery.h"

#define NANOSECTOSEC(elapsed_time) ((elapsed_time)/(double)1000000000)
#define BYTESTOMB(memory_cost) ((memory_cost)/(double)(1024 * 1024))
#define KBTOMB(memory_cost) ((memory_cost)/(double)(1024))

int main(int argc, char** argv) {
    MatchingCommand command(argc, argv);
    std::string input_query_graph_file = command.getQueryGraphFilePath();
    std::string input_data_graph_file = command.getDataGraphFilePath();
    std::string input_filter_type = command.getFilterType();
    std::string input_order_type = command.getOrderType();
    std::string input_engine_type = command.getEngineType();
    std::string input_max_embedding_num = command.getMaximumEmbeddingNum();
    std::string input_time_limit = command.getTimeLimit();

    /**
     * Output the command line information.
     */
    std::cout << "Command Line:" << std::endl;
    std::cout << "\tData Graph: " << input_data_graph_file << std::endl;
    std::cout << "\tQuery Graph: " << input_query_graph_file << std::endl;
    std::cout << "\tFilter Type: " << input_filter_type << std::endl;
    std::cout << "\tOrder Type: " << input_order_type << std::endl;
    std::cout << "\tEngine Type: " << input_engine_type << std::endl;
    std::cout << "\tOutput Limit: " << input_max_embedding_num << std::endl;
    std::cout << "\tTime Limit (seconds): " << input_time_limit << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;

    /**
     * Load input graphs.
     */
    std::cout << "Load graphs..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    Graph* query_graph = new Graph(true);
    query_graph->loadGraphFromFile(input_query_graph_file);
    query_graph->g_name = input_query_graph_file;
    query_graph->buildCoreTable();

    Graph* data_graph = new Graph(true);
    data_graph->loadGraphFromFile(input_data_graph_file);

    auto end = std::chrono::high_resolution_clock::now();

    double load_graphs_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    std::cout << "-----" << std::endl;
    std::cout << "Query Graph Meta Information" << std::endl;
    query_graph->printGraphMetaData();
    std::cout << "-----" << std::endl;
    data_graph->printGraphMetaData();

    std::cout << "--------------------------------------------------------------------" << std::endl;

    /**
     * Start queries.
     */

    std::cout << "Start queries..." << std::endl;
    std::cout << "-----" << std::endl;
    std::cout << "Filter candidates..." << std::endl;

    start = std::chrono::high_resolution_clock::now();

    ui** candidates = nullptr;
    ui* candidates_count = nullptr;
    ui* cfl_order = nullptr;
    TreeNode* cfl_tree = nullptr;
    ui* dpiso_order = nullptr;
    TreeNode* dpiso_tree = nullptr;
    std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> TE_Candidates;
    std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;
    if (input_filter_type == "LDF") {
        FilterVertices::LDFFilter(data_graph, query_graph, candidates, candidates_count);
    } else if (input_filter_type == "NLF") {
        FilterVertices::NLFFilter(data_graph, query_graph, candidates, candidates_count);
    } else if (input_filter_type == "CFL") {
        FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, cfl_order, cfl_tree);
    } else if (input_filter_type == "DPiso") {
        FilterVertices::DPisoFilter(data_graph, query_graph, candidates, candidates_count, dpiso_order, dpiso_tree);
    } else {
        std::cout << "The specified filter type '" << input_filter_type << "' is not supported." << std::endl;
        exit(-1);
    }

    // Sort the candidates to support the set intersections
    FilterVertices::sortCandidates(candidates, candidates_count, query_graph->getVerticesCount());

    end = std::chrono::high_resolution_clock::now();
    double filter_vertices_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    std::cout << "-----" << std::endl;
    std::cout << "Build indices..." << std::endl;

    start = std::chrono::high_resolution_clock::now();

    Edges ***edge_matrix = nullptr;
    edge_matrix = new Edges **[query_graph->getVerticesCount()];
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        edge_matrix[i] = new Edges *[query_graph->getVerticesCount()];
    }

    if (input_engine_type != "BSX") {
        BuildEdgeIndex::buildCansIdxIndex(data_graph, query_graph, candidates, candidates_count, edge_matrix);
    } else {
        BuildEdgeIndex::buildCansIndex(data_graph, query_graph, candidates, candidates_count, edge_matrix);
    }

    end = std::chrono::high_resolution_clock::now();
    double build_table_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    std::cout << "-----" << std::endl;
    std::cout << "Generate a matching order..." << std::endl;

    start = std::chrono::high_resolution_clock::now();

    ui* matching_order = nullptr;
    ui* pivots = nullptr;

    if (input_order_type == "GQL") {
        GenerateQueryPlan::generateGQLQueryPlan(data_graph, query_graph, candidates_count, matching_order, pivots);
    } else {
        std::cout << "The specified order type '" << input_order_type << "' is not supported." << std::endl;
    }

    end = std::chrono::high_resolution_clock::now();
    double generate_query_plan_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    GenerateQueryPlan::checkQueryPlanCorrectness(query_graph, matching_order, pivots);
    GenerateQueryPlan::printSimplifiedQueryPlan(query_graph, matching_order);

    std::cout << "-----" << std::endl;
    std::cout << "Enumerate..." << std::endl;
    size_t output_limit = 0;
    mpz_t embedding_cnt;
    if (input_max_embedding_num == "MAX") {
        output_limit = (size_t)-1;  // -1 means no limits
    } else {
        sscanf(input_max_embedding_num.c_str(), "%zu", &output_limit);
    }

    size_t call_count = 0;
    size_t time_limit = 0;
    size_t valid_vtx_count = 0;
    sscanf(input_time_limit.c_str(), "%zu", &time_limit);

    start = std::chrono::high_resolution_clock::now();

    if (input_engine_type == "BS1") {
        EvaluateQuery::BS1Engine(data_graph, query_graph, edge_matrix, candidates,
                                                      candidates_count, matching_order, pivots, output_limit, call_count, embedding_cnt);
    } else if (input_engine_type == "BSX") {
        EvaluateQuery::BSXEngine(data_graph, query_graph, edge_matrix, candidates, candidates_count,
                matching_order, output_limit, call_count, embedding_cnt);
    } else {
        std::cout << "The specified engine type '" << input_engine_type << "' is not supported." << std::endl;
        exit(-1);
    }

    end = std::chrono::high_resolution_clock::now();
    double enumeration_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Release memories..." << std::endl;
    /**
     * Release the allocated memories.
     */
    delete[] candidates_count;
    delete[] cfl_order;
    delete[] cfl_tree;
    delete[] dpiso_order;
    delete[] dpiso_tree;
    delete[] matching_order;
    delete[] pivots;
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        delete[] candidates[i];
    }
    delete[] candidates;

    if (edge_matrix != nullptr) {
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            for (ui j = 0; j < query_graph->getVerticesCount(); ++j) {
                delete edge_matrix[i][j];
            }
            delete[] edge_matrix[i];
        }
        delete[] edge_matrix;
    }

    delete query_graph;
    delete data_graph;

    /**
     * End.
     */
    std::cout << "--------------------------------------------------------------------" << std::endl;
    double preprocessing_time_in_ns = filter_vertices_time_in_ns + build_table_time_in_ns + generate_query_plan_time_in_ns;
    double total_time_in_ns = preprocessing_time_in_ns + enumeration_time_in_ns;

    printf("Load graphs time (seconds): %.4lf\n", NANOSECTOSEC(load_graphs_time_in_ns));
    printf("Filter vertices time (seconds): %.4lf\n", NANOSECTOSEC(filter_vertices_time_in_ns));
    printf("Build table time (seconds): %.4lf\n", NANOSECTOSEC(build_table_time_in_ns));
    printf("Generate query plan time (seconds): %.4lf\n", NANOSECTOSEC(generate_query_plan_time_in_ns));
    printf("Enumerate time (seconds): %.4lf\n", NANOSECTOSEC(enumeration_time_in_ns));
    printf("Preprocessing time (seconds): %.4lf\n", NANOSECTOSEC(preprocessing_time_in_ns));
    printf("Total time (seconds): %.4lf\n", NANOSECTOSEC(total_time_in_ns));
    gmp_printf("#Embeddings: %Zd\n", embedding_cnt);
    printf("Call Count: %zu\n", call_count);
    printf("Per Call Count Time (nanoseconds): %.4lf\n", enumeration_time_in_ns / (call_count == 0 ? 1 : call_count));
#ifdef ANALYZE_MEMORY
    printf("Memory cost (MB): %.4lf\n", KBTOMB(mem::getValue()));
#endif
    std::cout << "End." << std::endl;

    /**
     * Set the output stream and record the command line information
     */
    std::fstream output;
    output.open("./bsx_output.csv", std::ios::out | std::ios::app);

    output << input_query_graph_file;
    output << "," << input_data_graph_file;
    output << "," << input_filter_type;
    output << "," << input_order_type;
    output << "," << input_engine_type;
    output << "," << NANOSECTOSEC(load_graphs_time_in_ns);
    output << "," << NANOSECTOSEC(filter_vertices_time_in_ns);
    output << "," << NANOSECTOSEC(build_table_time_in_ns);
    output << "," << NANOSECTOSEC(generate_query_plan_time_in_ns);
    output << "," << NANOSECTOSEC(enumeration_time_in_ns);
    output << "," << NANOSECTOSEC(preprocessing_time_in_ns);
    output << "," << NANOSECTOSEC(total_time_in_ns);
    char *cnt = mpz_get_str(NULL, 10, embedding_cnt);
    output << "," << cnt;
    output << "," << call_count;
#ifdef ANALYZE_MEMORY
    output << "," << mem::getValue();  // KB
#endif
    output << std::endl;
    
    output.close();
    free(cnt);

    return 0;
}