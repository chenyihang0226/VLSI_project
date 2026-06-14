#ifndef SOLUTION_H
#define SOLUTION_H

#include <string>
#include "Graph.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>

using namespace std;

class Solution{
    public:
        // 读取 ICCAD 2021 TopoPart 格式的 benchmark 文件
    // 格式：Line 1 = 节点数, Line 2 = 线网数, 之后每行一个线网（0 索引节点）
    // 末尾部分：固定节点信息（第 i 行的节点固定到 FPGA i-1）
    // 读取后内部统一使用 1 索引节点编号
    void read_benchmark(Graph &graph, string benchmark_name,
                        vector<vector<int>> *fixed_lists = nullptr);

        /**
         * 多层次超图划分主算法（支持多路划分 + 多线程并行）。
         *
         * 多路划分：
         *   支持任意 k >= 2 路划分。k=2 时行为与二分完全一致。
         *   每个节点被分配到 0..k-1 中的一个 partition，割代价为涉及
         *   至少两个不同 partition 的线网数。
         *
         * 平衡约束：
         *   每路权重在 [total/k * (1-r), total/k * (0.5+r)] 范围内。
         *   r = kBalanceRatio，即标准形式 (1±eps) * W/k。
         *
         * 并行策略：
         *   外层 kRestarts 次重启通过 OpenMP #pragma omp parallel for 并行。
         *
         * @param graph            助教框架的图数据结构
         * @param X                输出（仅 k=2）：划分到 X 侧的顶点编号
         * @param Y                输出（仅 k=2）：划分到 Y 侧的顶点编号
         * @param num_threads      并行线程数。0 表示自动检测
         * @param k                划分路数（>= 2）
         * @param node_partition   输出（k>2 时）：node_id→partition 映射，1-indexed
         */
        // 返回值：最佳割代价（cutsize）
        int my_partition_algorithm(Graph graph, set<int> &X, set<int> &Y,
                                   unsigned int num_threads = 0,
                                   int k = 2,
                                   vector<int> *node_partition = nullptr,
                                   const string &topo_file = "",
                                   const vector<int> *fixed_of = nullptr);
};

#endif
