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
    // 读取后内部统一使用 1 索引节点编号
    void read_benchmark(Graph &graph, string benchmark_name);

        /**
         * 多层次超图二分划分主算法（支持多线程并行）。
         *
         * 并行策略：
         *   kRestarts 次随机重启是整个算法中最"尴尬并行"（embarrassingly parallel）
         *   的部分——每次重启独立构建粗化层次、独立运行 FM、独立维护最佳分区，
         *   彼此之间没有任何数据依赖。
         *
         *   实现上将 kRestarts 次重启分批（每批最多 num_threads 个），
         *   通过 std::async(launch::async, ...) 分发到多个线程并行执行。
         *   每批结束后汇总结果（std::mutex 保护全局 best_part/best_cut），
         *   再启动下一批。
         *
         * @param graph        助教框架的图数据结构（划分前已读取完毕）
         * @param X            输出参数：划分到 X 侧的顶点编号集合
         * @param Y            输出参数：划分到 Y 侧的顶点编号集合
         * @param num_threads  并行线程数。0 表示自动检测硬件并发度
         */
        void my_partition_algorithm(Graph graph, set<int> &X, set<int> &Y,
                                    unsigned int num_threads = 0);
};

#endif