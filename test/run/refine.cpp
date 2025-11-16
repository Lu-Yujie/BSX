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

void constructCommands() {
    std::string exefile = "/home/yujielu/subgraph/bsx/code/reviews/refinement/build/bin/BS";
    std::string con_path = "/home/yujielu/subgraph/bsx/code/reviews/refinement/test/conf";
    std::string file_suffix = ".txt";
    int time_limit = 1000;
    std::string max_num = "MAX";
    int num = 1000;
    int parallel_num = MAX_CONCURRENT_COMMANDS;
    std::string refinement_ratio = "25";  // 25, 50, 75, 100

    std::vector<std::string> graphs = {"citeseer", "dblp", "HPRD", "human", "maayan-figeys", "twitch", "web-Stanford", "wordnet-words", "YeastS", "youtube", "gplus"};
    std::vector<int> lsizes = {15, 30, 45, 60};
    std::vector<int> sizes = {10, 20, 30, 40, 50};

    for (const auto &graph : graphs) {
        std::string dataspace = "/var/lib/docker/subgraph/dataset/slabel/nlabel/" + graph;
        std::string outputspace = "/var/lib/docker/subgraph/output/bsx/empty/refine/r" + refinement_ratio + "/" + graph;

        for (const auto &lsize : lsizes) {
            std::string data_graph = "/L" + std::to_string(lsize) + "/" + graph + "-" + std::to_string(lsize) + file_suffix;
            for (const auto &size : sizes) {
                for (int j1 = 0; j1 <= num / parallel_num; ++j1) {
                    for (int j2 = 1; j2 <= parallel_num && j1 * parallel_num + j2 <= num; ++j2) {
                        int i = j1 * parallel_num + j2;
                        std::string query_graph = "/L" + std::to_string(lsize) + "/Q" + std::to_string(size) + "/Q" + std::to_string(size) + "-" + std::to_string(i) + file_suffix;
                        std::string output_file = "/L" + std::to_string(lsize) + "/Q" + std::to_string(size) + "/" + std::to_string(j2) + ".csv";

                        std::ostringstream command;
                        command << "timeout 30 " << exefile << " -d " << dataspace + data_graph << " -q " << dataspace + query_graph
                                << " -num " << max_num << " -time_limit " << time_limit << " -o " << outputspace + output_file
                                << " -conf " << con_path;

                        commands.push(command.str());
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
// ps -ef | grep -E "BS|refine" | grep -v grep | cut -c 9-16 | xargs kill -9

// clear old outputs
// find . -type f -name "*.csv" -delete
