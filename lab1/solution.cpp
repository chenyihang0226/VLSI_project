#include "solution.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <unordered_map>
#include <vector>

namespace {
constexpr double kBalanceRatio = 0.02;
constexpr int kMaxGreedyIters = 5;
// constexpr int kMaxSwapIters = 8;
// constexpr int kSwapSampleLimit = 200;

vector<int> collect_node_ids(Graph &graph) {
    vector<int> ids;
    for (const auto node : graph.get_nodes()) {
        ids.push_back(node->get_index());
    }
    return ids;
}

unordered_map<int, Node *> build_node_map(Graph &graph) {
    unordered_map<int, Node *> map;
    for (const auto node : graph.get_nodes()) {
        map[node->get_index()] = node;
    }
    return map;
}

void compute_balance_bounds(int n, int &min_x, int &max_x) {
    int min_val = static_cast<int>(ceil((0.5 - kBalanceRatio) * n));
    int max_val = static_cast<int>(floor((0.5 + kBalanceRatio) * n));
    if (min_val < 0) min_val = 0;
    if (max_val > n) max_val = n;
    min_x = min_val;
    max_x = max_val;
}

void init_random_balanced(Graph &graph, set<int> &X, set<int> &Y, mt19937 &rng) {
    X.clear();
    Y.clear();

    vector<int> ids = collect_node_ids(graph);
    shuffle(ids.begin(), ids.end(), rng);

    int n = static_cast<int>(ids.size());
    int min_x = 0, max_x = 0;
    compute_balance_bounds(n, min_x, max_x);
    int target_x = n / 2;
    if (target_x < min_x) target_x = min_x;
    if (target_x > max_x) target_x = max_x;

    for (int i = 0; i < n; ++i) {
        if (i < target_x) {
            X.insert(ids[i]);
        } else {
            Y.insert(ids[i]);
        }
    }
}

void build_net_counts(Graph &graph,
                      const set<int> &X,
                      vector<int> &count_x,
                      vector<int> &count_y) {
    int net_num = graph.get_net_num();
    count_x.assign(net_num, 0);
    count_y.assign(net_num, 0);

    for (const auto net : graph.get_nets()) {
        int idx = net->get_index();
        for (const auto node : net->get_nodes()) {
            int id = node->get_index();
            if (X.find(id) != X.end()) {
                count_x[idx]++;
            } else {
                count_y[idx]++;
            }
        }
    }
}

int delta_cut_for_move(const Node *node,
                       bool from_x,
                       const vector<int> &count_x,
                       const vector<int> &count_y) {
    int delta = 0;
    for (const auto net : node->get_nets()) {
        int idx = net->get_index();
        int cx = count_x[idx];
        int cy = count_y[idx];

        bool before_cut = (cx > 0 && cy > 0); //>0表示切开了
        if (from_x) {
            cx -= 1;
            cy += 1;
        } else {
            cy -= 1;
            cx += 1;
        }
        bool after_cut = (cx > 0 && cy > 0);
        delta += (after_cut ? 1 : 0) - (before_cut ? 1 : 0);
    }
    return delta;
}

void apply_move_counts(const Node *node,
                       bool from_x,
                       vector<int> &count_x,
                       vector<int> &count_y) {
    for (const auto net : node->get_nets()) {
        int idx = net->get_index();
        if (from_x) {
            count_x[idx] -= 1;
            count_y[idx] += 1;
        } else {
            count_y[idx] -= 1;
            count_x[idx] += 1;
        }
    }
}
}

void Solution::read_benchmark(Graph &graph, string benchmark_name) {
    ifstream file(benchmark_name);

    if(!file.is_open()) {
        cerr << "Failed to open the file!" << endl;
        exit(-1);
    }

    int edge_num, node_num;
    string line;
    getline(file >> ws, line);
    istringstream iss(line);
    iss >> edge_num;
    iss >> node_num;

    
    for(int i = 0; i < edge_num; i++) {
        getline(file, line);
        istringstream iss(line);
        int node_id;
        
        Net *net = graph.add_net(i);

        while(iss >> node_id) {
            Node *node = graph.get_or_create_node(node_id);
            node->add_net(net);
            net->add_node(node);
        }
        
    }

    file.close();
}

// Algorithm 1: Random balanced partition (baseline)
// void Solution::my_partition_algorithm(Graph graph, set<int> &X, set<int> &Y) {
//     mt19937 rng(42);
//     init_random_balanced(graph, X, Y, rng);
// }

// Algorithm 2: Greedy single-move improvement (default enabled)
void Solution::my_partition_algorithm(Graph graph, set<int> &X, set<int> &Y) {
    mt19937 rng(42);
    init_random_balanced(graph, X, Y, rng);

    int n = graph.get_node_num();
    int min_x = 0, max_x = 0;
    compute_balance_bounds(n, min_x, max_x);

    vector<int> count_x;
    vector<int> count_y;
    build_net_counts(graph, X, count_x, count_y);

    unordered_map<int, Node *> node_map = build_node_map(graph);

    for (int iter = 0; iter < kMaxGreedyIters; ++iter) {
        int best_delta = 0;
        int best_node = -1;
        bool best_from_x = true;

        int size_x = static_cast<int>(X.size());
        // int size_y = static_cast<int>(Y.size());

        for (const auto &pair : node_map) {
            int id = pair.first;
            const Node *node = pair.second;
            bool from_x = (X.find(id) != X.end());

            if (from_x) {
                if (size_x - 1 < min_x || size_x - 1 > max_x) {
                    continue;
                }
            } else {
                if (size_x + 1 < min_x || size_x + 1 > max_x) {
                    continue;
                }
            }

            int delta = delta_cut_for_move(node, from_x, count_x, count_y);
            if (delta < best_delta) {
                best_delta = delta;
                best_node = id;
                best_from_x = from_x;
            }
        }

        if (best_node == -1) {
            break;
        }

        Node *node = node_map[best_node];
        apply_move_counts(node, best_from_x, count_x, count_y);
        if (best_from_x) {
            X.erase(best_node);
            Y.insert(best_node);
        } else {
            Y.erase(best_node);
            X.insert(best_node);
        }
    }
}

// Algorithm 3: KL-style swap (simple sampled search)
// void Solution::my_partition_algorithm(Graph graph, set<int> &X, set<int> &Y) {
//     mt19937 rng(42);
//     init_random_balanced(graph, X, Y, rng);
//
//     vector<int> count_x;
//     vector<int> count_y;
//     build_net_counts(graph, X, count_x, count_y);
//
//     unordered_map<int, Node *> node_map = build_node_map(graph);
//
//     for (int iter = 0; iter < kMaxSwapIters; ++iter) {
//         vector<int> nodes_x(X.begin(), X.end());
//         vector<int> nodes_y(Y.begin(), Y.end());
//         shuffle(nodes_x.begin(), nodes_x.end(), rng);
//         shuffle(nodes_y.begin(), nodes_y.end(), rng);
//
//         if (static_cast<int>(nodes_x.size()) > kSwapSampleLimit) {
//             nodes_x.resize(kSwapSampleLimit);
//         }
//         if (static_cast<int>(nodes_y.size()) > kSwapSampleLimit) {
//             nodes_y.resize(kSwapSampleLimit);
//         }
//
//         int best_delta = 0;
//         int best_a = -1;
//         int best_b = -1;
//
//         for (int a : nodes_x) {
//             Node *node_a = node_map[a];
//             int delta_a = delta_cut_for_move(node_a, true, count_x, count_y);
//             apply_move_counts(node_a, true, count_x, count_y);
//
//             for (int b : nodes_y) {
//                 Node *node_b = node_map[b];
//                 int delta_b = delta_cut_for_move(node_b, false, count_x, count_y);
//                 int delta = delta_a + delta_b;
//                 if (delta < best_delta) {
//                     best_delta = delta;
//                     best_a = a;
//                     best_b = b;
//                 }
//             }
//
//             apply_move_counts(node_a, false, count_x, count_y);
//         }
//
//         if (best_a == -1) {
//             break;
//         }
//
//         Node *node_a = node_map[best_a];
//         Node *node_b = node_map[best_b];
//         apply_move_counts(node_a, true, count_x, count_y);
//         apply_move_counts(node_b, false, count_x, count_y);
//
//         X.erase(best_a);
//         Y.insert(best_a);
//         Y.erase(best_b);
//         X.insert(best_b);
//     }
// }