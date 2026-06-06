#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include "Net.h"
#include "Node.h"
#include "Graph.h"
#include "evaluate.h"
#include "solution.h"

using namespace std;

int main(int argc, char **argv) {

    Solution solution;

    if(argc != 2) {
        cout << "Usage: ./main benchmark_file" << endl;
        exit(-1);
    }
    string benchmark_name = argv[1];
    Graph graph;

    solution.read_benchmark(graph, benchmark_name);    

    cout << "Num nodes: " << graph.get_node_num() << endl;
    cout << "Num nets: " << graph.get_net_num() << endl;

    // TODO: 
    // 1. finish your partition algorithm
    // 2. output your partition result to a file
    // 3. evaluate your partition result

    // 05-15 solution

    set<int> X;
    set<int> Y;
    solution.my_partition_algorithm(graph, X, Y);

    // 从输入路径上生成分区文件
    string base_name = benchmark_name;
    size_t slash_pos = base_name.find_last_of('/');
    if (slash_pos != string::npos) {
        base_name = base_name.substr(slash_pos + 1);
    }
    size_t dot_pos = base_name.find_last_of('.');
    if (dot_pos != string::npos) {
        base_name = base_name.substr(0, dot_pos);
    }
    string partition_file = base_name + "_partition.txt";

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

    cout << evaluate(graph, partition_file) << endl;

    return 0;
}