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
  std::string exefile = "/home/yujielu/subgraph/bsx/code/BSX/build/bin/BS";
  std::string con_path = "/home/yujielu/subgraph/bsx/code/BSX/test/conf";
  std::string queryspace = "/var/lib/docker/subgraph/dataset/slabel/global";
  std::string file_suffix = ".txt";
  int time_limit = 300000;
  std::string max_num = "100000";

  std::vector<std::string> graphs = {"citeseer", "dblp", "HPRD", "human", "maayan-figeys", "twitch", "web-Stanford", "wordnet-words", "YeastS", "youtube"};
  std::vector<int> lsizes = {15, 30, 45, 60};
  std::vector<int> query_nums = {1, 2, 3, 4, 5, 6, 7, 8};  // different query graphs

  for (const auto &graph : graphs) {
    std::string dataspace = "/var/lib/docker/subgraph/dataset/slabel/nlabel/" + graph;
    std::string outputspace = "/var/lib/docker/subgraph/output/bsx/empty/global/o_10_5/" + graph;
    for (const auto &lsize : lsizes) {
      std::string data_graph = "/L" + std::to_string(lsize) + "/" + graph + "-" + std::to_string(lsize) + file_suffix;
      for (const auto &query_num : query_nums) {
        std::string query_graph = "/q" + std::to_string(query_num) + "/q" + std::to_string(query_num) + "-" + std::to_string(lsize) + file_suffix;
        std::string output_file = "/L" + std::to_string(lsize) + "/q" + std::to_string(query_num) + ".csv";

        std::ostringstream command;
        command << "timeout 2000 " << exefile << " -d " << dataspace + data_graph << " -q " << queryspace + query_graph
                << " -num " << max_num << " -time_limit " << time_limit << " -o " << outputspace + output_file
                << " -conf " << con_path;

        commands.push(command.str());
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
	constructCommands();
	std::vector<std::thread> threads;
	for (int i = 0; i < MAX_CONCURRENT_COMMANDS; ++i) {
		threads.emplace_back(workerThread);
	}

	for (auto &thread : threads) {
		thread.join();
	}

	std::cout << "All commands have been executed.\n";
	return 0;
}

// how to stop
// ps -ef | grep -E "BS|global" | grep -v grep | cut -c 9-16 | xargs kill -9

// clear old outputs
// find . -type f -name "*.csv" -delete
