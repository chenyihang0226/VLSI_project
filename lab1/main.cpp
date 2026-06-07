#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>
#include <string>
#include <set>
#include "Net.h"
#include "Node.h"
#include "Graph.h"
#include "evaluate.h"
#include "solution.h"

using namespace std;
namespace fs = std::filesystem;

int main(int argc, char **argv) {

    Solution solution;

    // ─── 命令行解析 ───────────────────────────────────────────────────────
    // 用法: ./main benchmark_file [num_threads]
    //   benchmark_file  —— 必须，指向 ICCAD 2021 case 文件的路径
    //   num_threads     —— 可选，并行线程数（默认 0 = 自动检测硬件并发度）
    // 示例:
    //   ./main "./ICCAD2021-TopoPart-Benchmarks/Generated Benchmarks/case1" 8
    if(argc < 2 || argc > 3) {
        cout << "Usage: ./main benchmark_file [num_threads]" << endl;
        exit(-1);
    }
    string benchmark_name = argv[1];

    unsigned int num_threads = 0; // 0 表示自动检测
    if (argc == 3) {
        num_threads = static_cast<unsigned int>(stoul(argv[2]));
    }

    Graph graph;

    // 读取 ICCAD 2021 格式的 benchmark 文件到 Graph（单线程 I/O 操作，无需并行）
    cout << "Reading benchmark file... " << flush;
    solution.read_benchmark(graph, benchmark_name);
    cout << "done." << endl;

    cout << "Num nodes: " << graph.get_node_num() << endl;
    cout << "Num nets: " << graph.get_net_num() << endl;

    // 若未指定线程数，自动检测 CPU 逻辑核心数
    if (num_threads == 0) {
        num_threads = thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 1;
    }
    cout << "Using " << num_threads << " thread(s)" << endl;

    // ─── 执行划分算法 ──────────────────────────────────────────────────
    // my_partition_algorithm 内部使用 std::async 将 kRestarts 次随机重启
    // 分发到 num_threads 个线程并行运行，或在 deadline (90s) 到达时停止.
    // 对于 case1 ~ case8 等大规模图，如时间预算不足可修改
    // kTotalTimeBudgetSeconds 重新编译。
    cout << "Running partition algorithm... " << flush;
    set<int> X;
    set<int> Y;
    solution.my_partition_algorithm(graph, X, Y, num_threads);
    cout << "done." << endl;

    // ─── 输出分区结果 ──────────────────────────────────────────────────
    // 从输入路径提取基准名，生成分区文件到 results/ 目录
    string base_name = benchmark_name;
    size_t slash_pos = base_name.find_last_of("/\\"); // 同时支持 / 和 \ 路径分隔符
    if (slash_pos != string::npos) {
        base_name = base_name.substr(slash_pos + 1);
    }
    size_t dot_pos = base_name.find_last_of('.');
    if (dot_pos != string::npos) {
        base_name = base_name.substr(0, dot_pos);
    }

    // 创建 results 目录（若不存在）
    fs::create_directory("results");
    string partition_file = "results/" + base_name + "_partition.txt";

    ofstream out(partition_file);
    if (!out.is_open()) {
        cerr << "Failed to open output file: " << partition_file << endl;
        return -1;
    }
    for (int i = 1; i <= graph.get_node_num(); ++i) {
        if (X.find(i) != X.end()) {
            out << 0 << "\n";
        } else {
            out << 1 << "\n";
        }
    }
    out.close();

    cout << "Output written to " << partition_file << endl;
    cout << "Cutsize: " << evaluate(graph, partition_file) << endl;

    return 0;
}