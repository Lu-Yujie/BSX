#include "Experiment.h"
#include "timeOp.h"

using namespace std;

#define NANOSECTOSEC(elapsed_time) ((elapsed_time)/(double)(1000*1000*1000))
#define BYTESTOMB(memory_cost) ((memory_cost)/(double)(1024 * 1024))

int main(int argc, char** argv) {
    MatchingCommand command(argc, argv);
    string input_query_graph_file = command.getQueryGraphFilePath();
    string input_data_graph_file = command.getDataGraphFilePath();
    string input_max_embedding_num = command.getMaximumEmbeddingNum();
    string input_time_limit = command.getTimeLimit();
    string input_output_file = command.getOutputFile();
    string input_conf_path = command.getConfPath();

    /**
     * Output the command line information.
     */
    std::cout << "Command Line:" << endl;
    std::cout << "\tData Graph: " << input_data_graph_file << endl;
    std::cout << "\tQuery Graph: " << input_query_graph_file << endl;
    std::cout << "\tOutput Limit: " << input_max_embedding_num << endl;
    std::cout << "\tTime Limit (milliseconds): " << input_time_limit << endl;
    std::cout << "\tOutput File Path: " << input_output_file << endl;
    std::cout << "\tConfiguration File Path: " << input_conf_path << endl;
    std::cout << "--------------------------------------------------------------------" << endl;

    /**
     * read methods info from configuration files
     */
    fstream output;
    output.open(input_output_file, ios::out | ios::app);
    Methods methods(input_conf_path);

    /**
     * Load input graphs.
     */
    std::cout << "Load graphs..." << endl;
    auto start = chrono::high_resolution_clock::now();

    Graph* query_graph = new Graph(true);
    query_graph->loadGraphFromFile(input_query_graph_file);
    query_graph->buildCoreTable();
    Graph* data_graph = new Graph(true);
    data_graph->loadGraphFromFile(input_data_graph_file);

    auto end = chrono::high_resolution_clock::now();
    auto load_graphs_time_in_ns = TimeOp::diffNan(start, end);

    std::cout << "-----" << endl;
    std::cout << "Query Graph Meta Information" << endl;
    query_graph->printGraphMetaData();
    std::cout << "-----" << endl;
    data_graph->printGraphMetaData();
    std::cout << "Query Graph Meta Information" << endl;
    std::cout << "--------------------------------------------------------------------" << endl;
    // write info to output file
    output << "Data Graph:" << input_data_graph_file << endl;
    output << "Query Graph:" << input_query_graph_file << endl;
    output << "load graph time:" << NANOSECTOSEC(load_graphs_time_in_ns) << endl;

    /**
     * init variables, set limits
     */
    // time
    auto lastSlash = input_query_graph_file.find_last_of('/');
    string query_name;
    if (lastSlash == std::string::npos) {
        query_name = input_query_graph_file;
    } else {
        query_name = input_query_graph_file.substr(lastSlash + 1);
    }
    int64_t time_limit; // 300s by default
    sscanf(input_time_limit.c_str(), "%ld", &time_limit); // second
    cout << time_limit << endl;
    // embedding_cnt
    uint64_t output_limit = 0;
    if (input_max_embedding_num == "MAX") {
        output_limit = numeric_limits<uint64_t>::max();
    }
    else {
        sscanf(input_max_embedding_num.c_str(), "%zu", &output_limit);
    }

    /**
     * Start query, scan methods
     */
    string filter_type, order_type, engine_type, combanition;
    mpz_t embedding_cnt;
    mpz_init(embedding_cnt);
    char cnt_buff[1000];
    // init all src_method first, can use virtual function here
    for (auto src : methods.src_method) {
        if (src == "VEQ") {
            VEQ::init_veq(input_data_graph_file, time_limit);
        } else {
            std::cout << "Please check conf file or put your method info here" << endl;
            exit(-1);
        }
    }
    do {
    // get methods
    if (methods.src_com) {  // src methods
        engine_type = methods.getCurMethods(methods.S_IDX);
        bool overtime = false;
        int64_t total_time = 0;
        if (engine_type == "VEQ") {
            overtime = VEQ::src_veq(input_query_graph_file, embedding_cnt, total_time, output_limit);
        }
        if (mpz_cmp_ui(embedding_cnt, 0) != 0) {
            mpz_get_str(cnt_buff, 10, embedding_cnt);
        } else {
            cnt_buff[0] = '0';
            cnt_buff[1] = 0;
        }
        output << query_name << "," << overtime << "," << engine_type << "," << engine_type << "," << engine_type << "," \
            << NANOSECTOSEC(total_time) << "," << cnt_buff << endl;
        continue;
    } else {  // combanition methods
        filter_type = methods.getCurMethods(methods.F_IDX);
        order_type = methods.getCurMethods(methods.O_IDX);
        engine_type = methods.getCurMethods(methods.E_IDX);
    }

    /************************************ init variables *****************************************/
    int64_t filter_time_in_ns = 0;
    int64_t build_table_time_in_ns = 0;
    int64_t order_time_in_ns = 0;
    int64_t engine_time_in_ns = 0;
    auto end_time = TimeOp::getClockNan();
    end_time += time_limit * 1000 * 1000;
    bool overtime = false;
    uint64_t call_cnt = 0;
    mpz_set_ui(embedding_cnt, 0);

    // filter variables
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    ui* tso_order = NULL;
    TreeNode* tso_tree = NULL;
    ui* cfl_order = NULL;
    TreeNode* cfl_tree = NULL;
    ui* dpiso_order = NULL;
    TreeNode* dpiso_tree = NULL;
    TreeNode* veq_tree = NULL;
    ui* veq_order = NULL;
    TreeNode* ceci_tree = NULL;
    ui* ceci_order = NULL;
    catalog* storage = NULL;
    vector<unordered_map<VertexID, vector<VertexID >>> TE_Candidates;
    vector<vector<unordered_map<VertexID, vector<VertexID>>>> NTE_Candidates;
    Edges ***edge_matrix = NULL;  // graph index

    // order variables
    ui* matching_order = NULL;
    ui* pivots = NULL;
    ui** weight_array = NULL;

    // engine variables
    uint64_t valid_vtx_count = 0;
    /************************************ end variables ******************************************/

    
    std::cout << "Start query: filter: " << filter_type << ", order: " << order_type
              << ", engine: " << engine_type << endl;

    std::cout << "filter part" << endl;
    start = chrono::high_resolution_clock::now();
    if (filter_type == "LDF") {
        FilterVertices::LDFFilter(data_graph, query_graph, candidates, candidates_count, end_time);
    } else if (filter_type == "NLF") {
        FilterVertices::NLFFilter(data_graph, query_graph, candidates, candidates_count, end_time);
    } else if (filter_type == "GQL") {
        FilterVertices::GQLFilter(data_graph, query_graph, candidates, candidates_count, end_time);
    } else if (filter_type == "TSO") {
        FilterVertices::TSOFilter(data_graph, query_graph, candidates, candidates_count, tso_order,
                                  tso_tree, end_time);
    } else if (filter_type == "CFL") {
        FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, cfl_order,
                                  cfl_tree, end_time);
    } else if (filter_type == "DPiso") {
        FilterVertices::DPisoFilter(data_graph, query_graph, candidates, candidates_count, dpiso_order,
                                    dpiso_tree, end_time);
    } else if (filter_type == "VEQ") {
        FilterVertices::VEQFilter(data_graph, query_graph, candidates, candidates_count, veq_order,
                                  veq_tree, end_time);
    } else if (filter_type == "CECI") {
        FilterVertices::CECIFilter(data_graph, query_graph, candidates, candidates_count, ceci_order,
                                   ceci_tree, TE_Candidates, NTE_Candidates, end_time);
    } else if (filter_type == "RM") {
        FilterVertices::RMFilter(data_graph, query_graph, candidates, candidates_count, storage, end_time);
    } else if (filter_type == "CaLiG") {
        FilterVertices::CaLiGFilter(data_graph, query_graph, candidates, candidates_count, end_time);
    } else if (filter_type == "null") {
        ;  // do nothing
    } else {
        std::cout << "The specified filter type '" << filter_type << "' is not supported." << endl;
    }
    end = chrono::high_resolution_clock::now();
    filter_time_in_ns = TimeOp::diffNan(start, end);
    if (TimeOp::getClockNan() > end_time) {
        overtime = true;
        goto EXIT;
    }

    // build edge matrix (index)
    start = chrono::high_resolution_clock::now();
    if (filter_type != "CECI" && filter_type != "null") { // indices of CECI don't have no-tree edges
        FilterVertices::sortCandidates(candidates, candidates_count, query_graph->getVerticesCount());
        edge_matrix = new Edges **[query_graph->getVerticesCount()];
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            edge_matrix[i] = new Edges *[query_graph->getVerticesCount()];
        }
        if (engine_type != "BSX") {
            BuildEdgeIndex::buildCansIdxIndex(data_graph, query_graph, candidates, candidates_count, edge_matrix);
        } else {
            BuildEdgeIndex::buildCansIndex(data_graph, query_graph, candidates, candidates_count, edge_matrix);
        }
    }
    end = chrono::high_resolution_clock::now();
    build_table_time_in_ns = TimeOp::diffNan(start, end);

    std::cout << "------------" << endl;
    std::cout << "Generate a matching order..." << endl;
    start = chrono::high_resolution_clock::now();
    if (order_type == "QSI") {
        GenerateQueryPlan::generateQSIQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots);
    } else if (order_type == "GQL") {
        GenerateQueryPlan::generateGQLQueryPlan(data_graph, query_graph, candidates_count, matching_order, pivots);
    } else if (order_type == "TSO") {
        if (tso_tree == NULL) {
            GenerateFilteringPlan::generateTSOFilterPlan(data_graph, query_graph, tso_tree, tso_order);
        }
        GenerateQueryPlan::generateTSOQueryPlan(query_graph, edge_matrix, matching_order, pivots, tso_tree, tso_order);
    } else if (order_type == "CFL") {
        if (cfl_tree == NULL) {
            int level_count;
            ui* level_offset;
            GenerateFilteringPlan::generateCFLFilterPlan(data_graph, query_graph, cfl_tree, cfl_order, level_count, level_offset);
            delete[] level_offset;
        }
        GenerateQueryPlan::generateCFLQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots, cfl_tree, cfl_order, candidates_count);
    } else if (order_type == "DPiso") {
        if (dpiso_tree == NULL) {
            GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, dpiso_tree, dpiso_order);
        }

        GenerateQueryPlan::generateDSPisoQueryPlan(query_graph, edge_matrix, matching_order, pivots, dpiso_tree, dpiso_order,
                                                    candidates_count, weight_array);
    } else if (order_type == "CECI") {
        GenerateQueryPlan::generateCECIQueryPlan(query_graph, ceci_tree, ceci_order, matching_order, pivots);
    } else if (order_type == "RI") {
        GenerateQueryPlan::generateRIQueryPlan(data_graph, query_graph, matching_order, pivots);
    } else if (order_type == "VF2PP") {
        GenerateQueryPlan::generateVF2PPQueryPlan(data_graph, query_graph, matching_order, pivots);
    } else if (order_type == "VF3") {
        GenerateQueryPlan::generateVF3QueryPlan(data_graph, query_graph, matching_order, pivots);
    } else if (order_type == "RM") {
        GenerateQueryPlan::generateRMQueryPlan(query_graph, matching_order, edge_matrix, pivots);
    } else if (order_type == "null") {
        ;  // do nothing
    } else {
        std::cout << "The specified order type '" << order_type << "' is not supported." << endl;
    }

    end = chrono::high_resolution_clock::now();
    order_time_in_ns = TimeOp::diffNan(start, end);

    if (order_type != "null") {
        GenerateQueryPlan::checkQueryPlanCorrectness(query_graph, matching_order, pivots);
        GenerateQueryPlan::printSimplifiedQueryPlan(query_graph, matching_order);
    }

    std::cout << "------------" << endl;
    std::cout << "Enumerate..." << endl;
    start = chrono::high_resolution_clock::now();
    if (engine_type == "EXPLORE") {
        overtime = EvaluateQuery::ExploreEngine(data_graph, query_graph, edge_matrix, candidates, candidates_count,
                                               matching_order, pivots, output_limit, call_cnt, embedding_cnt, end_time);
    } else if (engine_type == "LFTJ") {
        overtime = EvaluateQuery::LFTJ(data_graph, query_graph, edge_matrix, candidates, candidates_count,
                                       matching_order, output_limit, call_cnt, valid_vtx_count, embedding_cnt, end_time);
    } else if (engine_type == "GQL") {
        overtime = EvaluateQuery::GQLEngine(data_graph, query_graph, candidates, candidates_count,matching_order,
                                                      output_limit, call_cnt, embedding_cnt, end_time);
    } else if (engine_type == "QSI") {
        overtime = EvaluateQuery::QSIEngine(data_graph, query_graph, candidates, candidates_count, matching_order,
                                                      pivots, output_limit, call_cnt, embedding_cnt, end_time);
    } else if (engine_type == "VF3") {
        overtime = EvaluateQuery::VF3Engine(data_graph, query_graph, matching_order, pivots, output_limit, call_cnt,
                                                  embedding_cnt, end_time);
    } else if (engine_type == "VEQ") {
        if (veq_tree == NULL) {
            GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, veq_tree, veq_order);
        }
        overtime = EvaluateQuery::VEQEngine(data_graph, query_graph, veq_tree, edge_matrix, candidates, candidates_count,
                                                  output_limit, call_cnt, embedding_cnt, end_time);
    } else if (engine_type == "DPiso") {
        overtime = EvaluateQuery::DPisoEngine(data_graph, query_graph, dpiso_tree, edge_matrix, candidates,
                                                    candidates_count, weight_array, dpiso_order, output_limit, call_cnt,
                                                    embedding_cnt, end_time);
    } else if (engine_type == "CECI") {
        overtime = EvaluateQuery::CECIEngine(data_graph, query_graph, ceci_tree, candidates, candidates_count, TE_Candidates,
                                                   NTE_Candidates, ceci_order, output_limit, call_cnt, embedding_cnt, end_time);
    } else if (engine_type == "RM") {
        overtime = EvaluateQuery::RMEngine(query_graph, data_graph, storage, edge_matrix, candidates, candidates_count,
                                                        matching_order, output_limit, call_cnt, embedding_cnt, end_time);
    } else if (engine_type == "KSS") {
        overtime = EvaluateQuery::KSSEngine(query_graph, data_graph, edge_matrix, candidates, candidates_count,
                                                         matching_order, output_limit, call_cnt, embedding_cnt, end_time);
    } else if (engine_type == "BSX") {
        overtime = EvaluateQuery::BSXEngine(data_graph, query_graph, edge_matrix, candidates, candidates_count,
                                             output_limit, call_cnt, embedding_cnt, end_time);
    } else {
        std::cout << "The specified engine type '" << engine_type << "' is not supported." << endl;
    }

    end = chrono::high_resolution_clock::now();
    engine_time_in_ns = TimeOp::diffNan(start, end);

    EXIT:
    std::cout << "-----" << endl;
    std::cout << "Release memories..." << endl;
    delete[] candidates_count;
    delete[] tso_order;
    delete[] tso_tree;
    delete[] cfl_order;
    delete[] cfl_tree;
    delete[] dpiso_order;
    delete[] dpiso_tree;
    delete[] ceci_order;
    delete[] ceci_tree;
    delete[] matching_order;
    delete[] pivots;
    delete storage;
    if (candidates != NULL) {
        for (ui i = 0; i < query_graph->getVerticesCount(); i++) {
            delete[] candidates[i];
        }
        delete[] candidates;
    }
    if (edge_matrix != NULL) {
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            for (ui j = 0; j < query_graph->getVerticesCount(); ++j) {
                delete edge_matrix[i][j];
            }
            delete[] edge_matrix[i];
        }
        delete[] edge_matrix;
    }
    if (weight_array != NULL) {
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            delete[] weight_array[i];
        }
        delete[] weight_array;
    }

    // print info to file
    int64_t preprocessing_time_in_ns = filter_time_in_ns + build_table_time_in_ns + order_time_in_ns;
    int64_t total_time_in_ns = preprocessing_time_in_ns + engine_time_in_ns;
    if (mpz_cmp_ui(embedding_cnt, 0) != 0) {
        mpz_get_str(cnt_buff, 10, embedding_cnt);
    } else {
        cnt_buff[0] = '0';
        cnt_buff[1] = 0;
    }

    output << query_name << "," << overtime << "," << filter_type << "," << order_type << "," << engine_type << "," \
           << NANOSECTOSEC(total_time_in_ns) << "," << cnt_buff << endl;
    } while(methods.next());

    delete query_graph;
    delete data_graph;
    mpz_clear(embedding_cnt);
    output << "end" << endl;
    
    output.close();
    return 0;
}