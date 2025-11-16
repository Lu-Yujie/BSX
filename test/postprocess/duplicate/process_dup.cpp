#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <mutex>
#include <thread>

const int MAX_CONCURRENT_COMMANDS = 25;

std::queue<std::string> commands;
std::mutex queue_mutex;

void constructCommands() {
  std::string exefile = "/home/yujielu/subgraph/bsx/code/BSX/test/postprocess/duplicate/duplicate.py";
  std::string datapath = "/var/lib/docker/subgraph/output/bsx/duplicate/";
  int num = 100;
  std::vector<std::string> methods = {"bsx", "bs1"};
  std::vector<std::string> graphs = {"citeseer", "Figeys", "YeastS"};

  for (auto& graph : graphs) {
  for (auto& method : methods) {
  for (int i = 0; i < num; i++) {
      std::string command = "python " + exefile + " " + datapath + " " + graph + " " + method + " " + std::to_string(i);
      commands.push(command);
  }}}
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
// ps -ef | grep -E "duplicate.py|process_dup" | grep -v grep | cut -c 9-16 | xargs kill -9

// clear old outputs
// find . -type f -name "*.csv" -delete

