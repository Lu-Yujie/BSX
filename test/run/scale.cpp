#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <thread>
#include <queue>
#include <mutex>
#include <cstdlib>

const int MAX_CONCURRENT_COMMANDS = 40;

std::queue<std::string> commands;
std::mutex queue_mutex;

// Function to construct commands
void constructCommands() {
    std::string exefile = "/home/yujielu/subgraph/bsx/code/BSX/build/bin/BS";
    std::string con_path = "/home/yujielu/subgraph/bsx/code/BSX/test/conf";
    std::string max_num = "MAX";
    std::string file_suffix = ".txt";
    int time_limit = 1000;
    int num = 100;
    int parallel_num = MAX_CONCURRENT_COMMANDS;

    std::vector<std::string> graphs = {"human", "HPRD", "citeseer", "YeastS", "Figeys"};
    std::vector<std::string> scales = { "5", "10", "20", "50", "100"};
    std::vector<int> sizes = {10, 20, 30, 40, 50};

    for (const auto &graph : graphs) {
    for (const auto& scale : scales) {
        std::string dataspace = "/var/lib/docker/subgraph/dataset/slabel/nlabel/synthetic/EvoGraph/" + graph + "/" + graph + "_" + scale+ "/L30";
        std::string data_graph = dataspace + "/" + graph + "_" + scale + "-30" + file_suffix;
        std::string outputspace = "/var/lib/docker/subgraph/output/bsx/empty/EvoGraph/" + graph + "/" + graph + "_" + scale;

            for (const auto &size : sizes) {
                for (int j1 = 0; j1 <= num / parallel_num; ++j1) {
                    for (int j2 = 1; j2 <= parallel_num && j1 * parallel_num + j2 <= num; ++j2) {
                        int i = j1 * parallel_num + j2;
                        std::string query_graph = dataspace + "/Q" + std::to_string(size) + "/Q" + std::to_string(size) + "-" + std::to_string(i) + file_suffix;
                        std::string output_file = outputspace + "/Q" + std::to_string(size) + "/" + std::to_string(j2) + ".csv";

                        std::string command;
                        // command << "timeout 30 " << exepath << " -d " << data_graph << " -q " << query_graph
                        //                 << " -num " << max_num << " -time_limit " << time_limit << " -o " << output_file;
                        command = "timeout 20 " + exefile + " -d " + data_graph + " -q " + query_graph +
                                  " -time_limit " + std::to_string(time_limit) + " -o " + output_file +
                                  " -num " + max_num + " -conf " + con_path;
                        commands.push(command);
                    }
                }
            }
        }
    }
    return;
}

void workerThread() {
    while (true) {
        std::string command;
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            if (commands.empty()) {
                return;
            }
            command = std::move(commands.front());
            commands.pop();
        }
        int ret = std::system(command.c_str());
    }
}

int main() {
    // 将命令添加到队列
    constructCommands();

    // 创建并启动线程池
    std::vector<std::thread> threads;
    for (int i = 0; i < MAX_CONCURRENT_COMMANDS; ++i) {
        threads.emplace_back(workerThread);
    }

    // 等待所有线程结束
    for (auto &thread : threads) {
        thread.join();
    }

    std::cout << "All commands have been executed.\n";
    return 0;
}

// how to stop
// ps -ef | grep -E "BS|scale" | grep -v grep | cut -c 9-16 | xargs kill -9

// clear old outputs
// find . -type f -name "*.csv" -delete
