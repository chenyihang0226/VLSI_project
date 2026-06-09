#include <chrono>
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
    // 用法: ./main benchmark_file [num_threads] [k] [--topo topo_file]
    // 示例:
    //   ./main ".../case1"                    # 二分
    //   ./main ".../case1" 8 8                # 8 路
    //   ./main ".../case1" 8 43 --topo "FPGA Graph/MFS2"  # 43 路 + 拓扑约束
    if(argc < 2) {
        cout << "Usage: ./main benchmark_file [num_threads] [k] [--topo topo_file]" << endl;
        exit(-1);
    }
    string benchmark_name = argv[1];

    unsigned int num_threads = 0;
    int k = 2;
    string topo_file = "";

    // Parse positional args and optional --topo flag
    int pos_args = 0;
    for (int i = 2; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--topo" && i + 1 < argc) {
            topo_file = argv[++i];
        } else if (pos_args == 0) {
            num_threads = static_cast<unsigned int>(stoul(arg));
            pos_args++;
        } else if (pos_args == 1) {
            k = stoi(arg);
            if (k < 2) { cerr << "k must be >= 2" << endl; return -1; }
            pos_args++;
        } else {
            cerr << "Unexpected argument: " << arg << endl;
            return -1;
        }
    }

    Graph graph;

    // 读取网表（若指定了拓扑约束，同时读取固定节点信息）
    vector<vector<int>> fixed_lists;
    vector<int> fixed_of;  // node_id → partition, -1 = unfixed
    bool has_fixed = false;

    cout << "Reading benchmark file... " << flush;
    if (!topo_file.empty()) {
        solution.read_benchmark(graph, benchmark_name, &fixed_lists);
    } else {
        solution.read_benchmark(graph, benchmark_name);
    }
    cout << "done." << endl;

    cout << "Num nodes: " << graph.get_node_num() << endl;
    cout << "Num nets: " << graph.get_net_num() << endl;
    cout << "k-way: " << k << endl;

    if (!topo_file.empty()) {
        cout << "Topology: " << topo_file << endl;
    }

    if (num_threads == 0) {
        num_threads = thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 1;
    }
    cout << "Using " << num_threads << " thread(s)" << endl;

    // 构建固定节点映射（fixed_of[node_id] = partition_id, 1-indexed）
    if (!fixed_lists.empty()) {
        fixed_of.assign(graph.get_node_num() + 1, -1);
        for (size_t p = 0; p < fixed_lists.size(); ++p) {
            for (int node_id : fixed_lists[p]) {
                if (node_id < (int)fixed_of.size())
                    fixed_of[node_id] = static_cast<int>(p);
            }
        }
        int fixed_count = 0;
        for (int f : fixed_of) if (f != -1) fixed_count++;
        cout << "Fixed nodes: " << fixed_count << " (across " << fixed_lists.size() << " FPGAs)" << endl;
        has_fixed = true;
    }

    // ─── 执行划分算法（计时） ────────────────────────────────────────────
    cout << "Running partition algorithm... " << flush;
    auto algo_start = chrono::steady_clock::now();
    set<int> X, Y;
    vector<int> node_partition; // node_id (1-indexed) → partition (0..k-1)
    const vector<int> *fixed_ptr = has_fixed ? &fixed_of : nullptr;
    int best_cut = solution.my_partition_algorithm(graph, X, Y, num_threads, k,
                                                    &node_partition,
                                                    topo_file, fixed_ptr);
    auto algo_end = chrono::steady_clock::now();
    double algo_seconds = chrono::duration<double>(algo_end - algo_start).count();
    cout << "done. (" << algo_seconds << " s)" << endl;
    cout << "Best cutsize found: " << best_cut << endl;

    // ─── 输出分区结果 ──────────────────────────────────────────────────
    string base_name = benchmark_name;
    size_t slash_pos = base_name.find_last_of("/\\");
    if (slash_pos != string::npos) {
        base_name = base_name.substr(slash_pos + 1);
    }
    size_t dot_pos = base_name.find_last_of('.');
    if (dot_pos != string::npos) {
        base_name = base_name.substr(0, dot_pos);
    }

    fs::create_directory("results");
    string partition_file = "results/" + base_name + "_partition.txt";

    ofstream out(partition_file);
    if (!out.is_open()) {
        cerr << "Failed to open output file: " << partition_file << endl;
        return -1;
    }

    if (k == 2) {
        // k=2: 输出 0/1，兼容旧的 evaluate
        for (int i = 1; i <= graph.get_node_num(); ++i) {
            out << (X.find(i) != X.end() ? 0 : 1) << "\n";
        }
    } else {
        // k>2: 输出 partition ID (0..k-1)
        for (int i = 1; i <= graph.get_node_num(); ++i) {
            out << node_partition[i] << "\n";
        }
    }
    out.close();

    cout << "Output written to " << partition_file << endl;
    cout << "Cutsize: " << best_cut << "  Time: " << algo_seconds << " s" << endl;

    return 0;
}