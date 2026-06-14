#include "solution.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <omp.h>      // OpenMP: #pragma omp parallel for, omp_get_max_threads()
#include <queue>
#include <random>
#include <unordered_map>
#include <vector>

namespace {
constexpr double kBalanceRatio = 0.02; //平衡约束参数
constexpr int kRestarts = 8; //针对启发式优化，不同重启产生不同的搜索路径
constexpr int kMaxLevels = 5;
constexpr int kMinCoarsestNodes = 320;
constexpr double kStopCoarsenRatio = 0.92;

// 局部最优更新参数：针对最细层以及不同粗化层之间设置不同的pass、step以及candidate
// 其中最细层参数最激进，最粗层参数最平缓
constexpr int kCoarseFmPasses = 8;
constexpr int kMidFmPasses = 8;
constexpr int kFineFmPasses = 10;
constexpr int kCoarseMaxSteps = 9000;
constexpr int kMidMaxSteps = 10000;
constexpr int kFineMaxSteps = 14000;
constexpr int kCoarseCandidateSize = 1800;
constexpr int kMidCandidateSize = 2200;
constexpr int kFineCandidateSize = 3000;
constexpr int kCoarsePolishPasses = 2;
constexpr int kMidPolishPasses = 2;
constexpr int kFinePolishPasses = 3;
constexpr int kCoarsePolishMaxSteps = 2500;
constexpr int kMidPolishMaxSteps = 3000;
constexpr int kFinePolishMaxSteps = 4500;
constexpr int kCoarsePolishCandidate = 2500;
constexpr int kMidPolishCandidate = 3200;
constexpr int kFinePolishCandidate = 4500;
constexpr int kTopoLambdaCoarse = 1;
constexpr int kTopoLambdaMid = 2;
constexpr int kTopoLambdaFine = 4;

// 时间预算（基准值，实际在 my_partition_algorithm 中按图规模自适应）
constexpr int kBudgetBaseTime = 90;        // 基准时间（秒），对应 ibm01 约 12k 节点
constexpr int kBudgetBaseNodes = 12000;    // 基准节点数
constexpr double kBudgetExponent = 0.8;    // 规模缩放指数（<1 为次线性）
constexpr int kBudgetMinSeconds = 60;      // 最小预算
constexpr int kBudgetMaxSeconds = 3600;    // 最大预算（1 小时）

// 粗化过程关键变量
constexpr int kHeavyNetLimit = 64; //用于限制粗化过程中net的选择
constexpr double kMatchDegreePenalty = 0.06; //控制候选节点评分，保证候选节点度数越高，评分越低
constexpr double kSmallNetBoost = 1.7; //如果一个net大小较小就加权，因为net越小越能体现强相关
constexpr int kInitModeCount = 3;

struct Hypergraph {
    int n = 0;
    vector<vector<int>> nets; // 连接 i 线网上的所有节点的内部索引
    vector<vector<int>> node_nets; // i 节点所连接的所有线网的索引
    vector<int> node_weight;
};

// 拓扑图结构：用于拓扑约束（Topology Constraints）
struct TopoGraph {
    int n = 0;                       // FPGA 节点数
    vector<bool> adj_matrix;         // 邻接矩阵 (n*n)，O(1) 查询连通性
};

// 读取拓扑图文件（格式：第一行 "节点数 边数"，后续每行 "u v"）
TopoGraph read_topo_graph(const string &filename) {
    ifstream file(filename);
    TopoGraph g;
    int m;
    if (!(file >> g.n >> m)) {
        cerr << "Failed to read topology file: " << filename << endl;
        exit(-1);
    }
    g.adj_matrix.assign(g.n * g.n, false);
    for (int i = 0; i < g.n; ++i)
        g.adj_matrix[i * g.n + i] = true;  // 自环
    for (int i = 0; i < m; ++i) {
        int u, v;
        if (!(file >> u >> v)) break;
        g.adj_matrix[u * g.n + v] = true;
        g.adj_matrix[v * g.n + u] = true;
    }
    return g;
}

// 加速版拓扑约束检查：利用 net 当前占用的 partition 集合，避免遍历所有 pin。
[[maybe_unused]] bool topo_check_move_fast(const TopoGraph &topo, const Hypergraph &hg,
                          int u, int to_part, const vector<int> &part,
                          const vector<int> &count,
                          const vector<int> &net_parts,
                          const vector<int> &net_parts_size,
                          int k) {
    int from_part = part[u];
    int n = topo.n;
    for (int net_idx : hg.node_nets[u]) {
        int base = net_idx * k;
        // to_part 已在该 net 中：增加同分区 pin 不会产生新的跨分区连接。
        if (count[base + to_part] > 0) continue;

        int sz = net_parts_size[net_idx];
        for (int i = 0; i < sz; ++i) {
            int p = net_parts[base + i];
            if (p == to_part) continue;
            // 移动后 from_part 会从该 net 消失，无需检查与其的邻接关系。
            if (p == from_part && count[base + from_part] == 1) continue;
            if (!topo.adj_matrix[to_part * n + p]) return false;
        }
    }
    return true;
}

// 拓扑约束检查（初始化中使用）：只检查已放置的邻居
[[maybe_unused]] bool topo_check_init(const TopoGraph &topo, const Hypergraph &hg,
                     int u, int to_part, const vector<int> &part) {
    for (int net_idx : hg.node_nets[u]) {
        for (int v : hg.nets[net_idx]) {
            if (v == u || part[v] == -1) continue;
            if (part[v] != to_part) {
                if (!topo.adj_matrix[to_part * topo.n + part[v]])
                    return false;
            }
        }
    }
    return true;
}

int topo_init_violation_count(const TopoGraph &topo, const Hypergraph &hg,
                              int u, int to_part, const vector<int> &part) {
    int violations = 0;
    for (int net_idx : hg.node_nets[u]) {
        bool bad_net = false;
        for (int v : hg.nets[net_idx]) {
            if (v == u || part[v] == -1) continue;
            int p = part[v];
            if (p != to_part && !topo.adj_matrix[to_part * topo.n + p]) {
                bad_net = true;
                break;
            }
        }
        if (bad_net) violations++;
    }
    return violations;
}

struct CoarsenResult {
    Hypergraph coarse;
    vector<int> fine_to_coarse;
    vector<int> coarse_fixed;
};

struct MoveRecord {
    int node = -1;
    int from_part = 0;   // 移动前的 partition (0..k-1)
    int to_part = 0;     // 移动目标 partition (0..k-1, != from_part)
    int gain = 0;
};

struct GainEntry {
    int gain = 0;
    int node = -1;
    int to_part = 0;     // 移动目标 partition
    int version = 0;
};

struct GainEntryCmp {
    bool operator()(const GainEntry &a, const GainEntry &b) const {
        return a.gain < b.gain;
    }
};

void compute_balance_bounds_kway(int total_weight, int k,
                                  vector<int> &min_w, vector<int> &max_w) {
    // k 路平衡约束：
    //   每路的目标权重 ≈ total_weight / k
    //   允许偏差为每份目标的 kBalanceRatio 倍，即标准形式 [(1-eps)*W/k, (1+eps)*W/k]
    //   这样 k=2 时与旧的 48%-52% 兼容，k=43 时不会允许某路几乎为空。
    double per_part = static_cast<double>(total_weight) / static_cast<double>(k);
    double margin = kBalanceRatio * per_part;
    int min_val = max(0, static_cast<int>(ceil(per_part - margin)));
    int max_val = min(total_weight, static_cast<int>(floor(per_part + margin)));
    min_w.assign(k, min_val);
    max_w.assign(k, max_val);
}

void rebuild_node_nets(Hypergraph &hg) {
    hg.node_nets.assign(hg.n, {});
    for (int net_idx = 0; net_idx < static_cast<int>(hg.nets.size()); ++net_idx) {
        for (int u : hg.nets[net_idx]) {
            hg.node_nets[u].push_back(net_idx);
        }
    }
}

Hypergraph build_fine_hypergraph(Graph &graph, vector<int> &idx_to_id) {
    Hypergraph hg;

    // 重新建立索引，编号为 0,1,2
    vector<int> ids;
    ids.reserve(graph.get_node_num());
    for (const auto node : graph.get_nodes()) {
        ids.push_back(node->get_index());
    }
    sort(ids.begin(), ids.end());
    idx_to_id = ids;

    unordered_map<int, int> id_to_idx; // 建立 "原始ID -> 内部索引" 的哈希表
    id_to_idx.reserve(ids.size() * 2 + 1);
    for (int i = 0; i < static_cast<int>(ids.size()); ++i) {
        id_to_idx[ids[i]] = i;
    }

    hg.n = static_cast<int>(ids.size()); // 节点数量
    hg.node_weight.assign(hg.n, 1); // 基础状态下每个节点的权重都是 1

    // 建立 线 - 点 网
    hg.nets.reserve(graph.get_net_num());
    for (const auto net : graph.get_nets()) {
        vector<int> nodes;
        nodes.reserve(net->get_nodes().size());
        for (const auto node : net->get_nodes()) {
            auto it = id_to_idx.find(node->get_index());
            if (it != id_to_idx.end()) {
                nodes.push_back(it->second);
            }
        }
        sort(nodes.begin(), nodes.end());
        nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end()); // 防止脏数据有重复

        // 过滤只有一个节点的无用线网
        if (nodes.size() >= 2) {
            hg.nets.push_back(std::move(nodes));
        }
    }

    // 构造 点 - 线 网
    rebuild_node_nets(hg);
    return hg;
}

// 粗化，V 的左边
bool fixed_compatible(const vector<int> *fixed_of, int a, int b) {
    if (!fixed_of) return true;
    int fa = (*fixed_of)[a];
    int fb = (*fixed_of)[b];
    return fa == -1 || fb == -1 || fa == fb;
}

CoarsenResult coarsen_once(const Hypergraph &fine, mt19937 &rng,
                           const vector<int> *fixed_of = nullptr) {
    CoarsenResult result;
    int n = fine.n;
    result.fine_to_coarse.assign(n, -1);

    vector<char> unmatched(n, 1);
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    // 随机种子打乱节点顺序
    shuffle(order.begin(), order.end(), rng);
    // 稳定排序，节点连接的线网多的排前面（度数高的对结果的影响更大）
    stable_sort(order.begin(), order.end(), [&fine](int a, int b) {
        return fine.node_nets[a].size() > fine.node_nets[b].size();
    });

    vector<vector<int>> groups;
    groups.reserve(n);

    for (int u : order) {
        if (!unmatched[u]) {
            continue;
        }
        unmatched[u] = 0;

        int best_v = -1;
        double best_score = -1.0;
        unordered_map<int, double> score;
        score.reserve(64);

        for (int net_idx : fine.node_nets[u]) {
            const auto &nodes = fine.nets[net_idx];
            int net_size = static_cast<int>(nodes.size());
            // 过滤掉连接了太多节点的线网，因为小net更能够说明节点之间关系紧密
            if (net_size > kHeavyNetLimit || net_size <= 1) {
                continue;
            }
            // 线网越小，权重越高
            double net_weight = 1.0 / static_cast<double>(net_size - 1);
            // 奖励
            if (net_size <= 4) {
                net_weight *= kSmallNetBoost;
            }
            for (int v : nodes) {
                if (v == u || !unmatched[v]) {
                    continue;
                }
                if (!fixed_compatible(fixed_of, u, v)) {
                    continue;
                }
                // 为了防止“超级枢纽“ 的出现，这里设置惩罚机制
                double deg_pen = 1.0 + kMatchDegreePenalty * static_cast<double>(fine.node_nets[v].size());
                // 得分：线网权重 / 惩罚系数
                double s = (score[v] += (net_weight / deg_pen));
                if (s > best_score) {
                    best_score = s;
                    best_v = v;
                } else if (best_v != -1 && fabs(s - best_score) < 1e-12) {
                    if (fine.node_nets[v].size() < fine.node_nets[best_v].size()) {
                        // 分数相同的，优先选择度数小的
                        best_v = v;
                    }
                }
            }
        }
        if (best_v != -1) {
            unmatched[best_v] = 0;
            // 好朋友组队
            groups.push_back({u, best_v});
        } else {
            groups.push_back({u});
        }
    }

    Hypergraph coarse;
    coarse.n = static_cast<int>(groups.size());
    coarse.node_weight.assign(coarse.n, 0);
    result.coarse_fixed.assign(coarse.n, -1);

    for (int cid = 0; cid < coarse.n; ++cid) {
        for (int u : groups[cid]) {
            result.fine_to_coarse[u] = cid; // 记录映射关系
            coarse.node_weight[cid] += fine.node_weight[u]; // 权重加和
        }
    }

    if (fixed_of) {
        for (int u = 0; u < n; ++u) {
            int p = (*fixed_of)[u];
            if (p != -1) {
                result.coarse_fixed[result.fine_to_coarse[u]] = p;
            }
        }
    }

    coarse.nets.reserve(fine.nets.size());
    for (const auto &net : fine.nets) {
        vector<int> coarse_nodes;
        coarse_nodes.reserve(net.size());
        for (int u : net) {
            coarse_nodes.push_back(result.fine_to_coarse[u]);
        }
        sort(coarse_nodes.begin(), coarse_nodes.end());
        // 合并后可能会有重复，都合并到同个 cid ，所以去重
        coarse_nodes.erase(unique(coarse_nodes.begin(), coarse_nodes.end()), coarse_nodes.end());
        // 如果一条线网的 size 变成了 1，说明只有一个节点了，就不进去作为网了
        if (coarse_nodes.size() >= 2) {
            coarse.nets.push_back(std::move(coarse_nodes));
        }
    }

    // 构造 点 - 线 网
    rebuild_node_nets(coarse);
    result.coarse = std::move(coarse);
    return result;
}

// ========== k-way initial partition methods ==========

// 1. Degree-balanced: sort by degree (descending), assign to lightest partition
vector<int> init_degree_balanced_kway(const Hypergraph &hg, int k, mt19937 &rng,
                                       const TopoGraph *topo,
                                       const vector<int> *fixed_of) {
    vector<int> part(hg.n, -1);
    vector<int> order(hg.n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);
    stable_sort(order.begin(), order.end(), [&hg](int a, int b) {
        return hg.node_nets[a].size() > hg.node_nets[b].size();
    });

    int total_weight = accumulate(hg.node_weight.begin(), hg.node_weight.end(), 0);
    vector<int> min_w(k), max_w(k);
    compute_balance_bounds_kway(total_weight, k, min_w, max_w);

    vector<int> weight(k, 0);

    // Pre-assign fixed nodes
    if (fixed_of) {
        for (int u = 0; u < hg.n; ++u) {
            int p = (*fixed_of)[u];
            if (p != -1) {
                if (p < 0 || p >= k) return {};
                part[u] = p;
                weight[p] += hg.node_weight[u];
                if (weight[p] > max_w[p]) return {};
            }
        }
    }

    for (int u : order) {
        if (part[u] != -1) continue;  // skip fixed nodes
        int w = hg.node_weight[u];
        int best_p = -1;
        int lightest = numeric_limits<int>::max();
        for (int p = 0; p < k; ++p) {
            if (weight[p] + w > max_w[p]) continue;
            if (topo && !topo_check_init(*topo, hg, u, p, part)) continue;
            if (weight[p] < lightest) { lightest = weight[p]; best_p = p; }
        }
        if (best_p == -1) return {};
        part[u] = best_p;
        weight[best_p] += w;
    }
    return part;
}

// 2. Random-balanced: random order, assign to lightest feasible partition
vector<int> init_random_balanced_kway(const Hypergraph &hg, int k, mt19937 &rng,
                                       const TopoGraph *topo,
                                       const vector<int> *fixed_of) {
    vector<int> part(hg.n, -1);
    vector<int> order(hg.n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);

    int total_weight = accumulate(hg.node_weight.begin(), hg.node_weight.end(), 0);
    vector<int> min_w(k), max_w(k);
    compute_balance_bounds_kway(total_weight, k, min_w, max_w);

    vector<int> weight(k, 0);

    if (fixed_of) {
        for (int u = 0; u < hg.n; ++u) {
            int p = (*fixed_of)[u];
            if (p != -1) {
                if (p < 0 || p >= k) return {};
                part[u] = p;
                weight[p] += hg.node_weight[u];
                if (weight[p] > max_w[p]) return {};
            }
        }
    }

    for (int u : order) {
        if (part[u] != -1) continue;
        int w = hg.node_weight[u];
        int best_p = -1;
        int lightest = numeric_limits<int>::max();
        for (int p = 0; p < k; ++p) {
            if (weight[p] + w > max_w[p]) continue;
            if (topo && !topo_check_init(*topo, hg, u, p, part)) continue;
            if (weight[p] < lightest) { lightest = weight[p]; best_p = p; }
        }
        if (best_p == -1) return {};
        part[u] = best_p;
        weight[best_p] += w;
    }
    return part;
}

// 3. Small-net-weighted: nodes in small nets get higher score, then assign
vector<int> init_smallnet_balanced_kway(const Hypergraph &hg, int k, mt19937 &rng,
                                         const TopoGraph *topo,
                                         const vector<int> *fixed_of) {
    vector<double> score(hg.n, 0.0);
    for (const auto &net : hg.nets) {
        int net_size = static_cast<int>(net.size());
        if (net_size <= 1) continue;
        double w = 1.0 / static_cast<double>(net_size - 1);
        if (net_size <= 4) w *= kSmallNetBoost;
        for (int u : net) score[u] += w;
    }

    vector<int> order(hg.n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);
    stable_sort(order.begin(), order.end(), [&score](int a, int b) {
        if (fabs(score[a] - score[b]) > 1e-12) return score[a] > score[b];
        return a < b;
    });

    vector<int> part(hg.n, -1);
    int total_weight = accumulate(hg.node_weight.begin(), hg.node_weight.end(), 0);
    vector<int> min_w(k), max_w(k);
    compute_balance_bounds_kway(total_weight, k, min_w, max_w);

    vector<int> weight(k, 0);
    vector<double> score_sum(k, 0.0);

    if (fixed_of) {
        for (int u = 0; u < hg.n; ++u) {
            int p = (*fixed_of)[u];
            if (p != -1) {
                if (p < 0 || p >= k) return {};
                part[u] = p;
                weight[p] += hg.node_weight[u];
                if (weight[p] > max_w[p]) return {};
            }
        }
    }

    for (int u : order) {
        if (part[u] != -1) continue;
        int w = hg.node_weight[u];
        int best_p = -1;
        int lightest = numeric_limits<int>::max();
        for (int p = 0; p < k; ++p) {
            if (weight[p] + w > max_w[p]) continue;
            if (topo && !topo_check_init(*topo, hg, u, p, part)) continue;
            if (weight[p] < lightest) { lightest = weight[p]; best_p = p; }
        }
        if (best_p == -1) {
            double min_score = numeric_limits<double>::max();
            for (int p = 0; p < k; ++p) {
                if (weight[p] + w <= max_w[p] && score_sum[p] < min_score) {
                    min_score = score_sum[p]; best_p = p;
                }
            }
            if (best_p == -1) return {};
        }
        part[u] = best_p;
        weight[best_p] += w;
        score_sum[best_p] += score[u];
    }
    return part;
}

vector<int> init_soft_topology_balanced_kway(const Hypergraph &hg, int k,
                                             mt19937 &rng,
                                             const TopoGraph *topo,
                                             const vector<int> *fixed_of,
                                             int mode) {
    if (!topo) {
        return init_random_balanced_kway(hg, k, rng, nullptr, fixed_of);
    }

    vector<double> smallnet_score(hg.n, 0.0);
    if (mode == 2) {
        for (const auto &net : hg.nets) {
            int net_size = static_cast<int>(net.size());
            if (net_size <= 1) continue;
            double w = 1.0 / static_cast<double>(net_size - 1);
            if (net_size <= 4) w *= kSmallNetBoost;
            for (int u : net) smallnet_score[u] += w;
        }
    }

    vector<int> order(hg.n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);
    if (mode == 0) {
        stable_sort(order.begin(), order.end(), [&hg](int a, int b) {
            return hg.node_nets[a].size() > hg.node_nets[b].size();
        });
    } else if (mode == 2) {
        stable_sort(order.begin(), order.end(), [&smallnet_score](int a, int b) {
            if (fabs(smallnet_score[a] - smallnet_score[b]) > 1e-12) {
                return smallnet_score[a] > smallnet_score[b];
            }
            return a < b;
        });
    }

    int total_weight = accumulate(hg.node_weight.begin(), hg.node_weight.end(), 0);
    vector<int> min_w(k), max_w(k);
    compute_balance_bounds_kway(total_weight, k, min_w, max_w);

    vector<int> part(hg.n, -1);
    vector<int> weight(k, 0);
    if (fixed_of) {
        for (int u = 0; u < hg.n; ++u) {
            int p = (*fixed_of)[u];
            if (p != -1) {
                if (p < 0 || p >= k) return {};
                part[u] = p;
                weight[p] += hg.node_weight[u];
                if (weight[p] > max_w[p]) return {};
            }
        }
    }

    auto cut_proxy = [&](int u, int to_part) {
        int proxy = 0;
        int affinity = 0;
        for (int net_idx : hg.node_nets[u]) {
            bool has_other = false;
            bool has_same = false;
            for (int v : hg.nets[net_idx]) {
                if (v == u || part[v] == -1) continue;
                if (part[v] == to_part) has_same = true;
                else has_other = true;
            }
            if (has_other) proxy++;
            if (has_same) affinity++;
        }
        return pair<int, int>(proxy, affinity);
    };

    for (int u : order) {
        if (part[u] != -1) continue;
        int w = hg.node_weight[u];
        int best_p = -1;
        int best_topo = numeric_limits<int>::max();
        int best_cut_proxy = numeric_limits<int>::max();
        int best_affinity = numeric_limits<int>::min();
        int best_weight = numeric_limits<int>::max();

        for (int p = 0; p < k; ++p) {
            if (weight[p] + w > max_w[p]) continue;
            int topo_bad = topo_init_violation_count(*topo, hg, u, p, part);
            auto cp = cut_proxy(u, p);
            int cproxy = cp.first;
            int affinity = cp.second;
            if (topo_bad < best_topo ||
                (topo_bad == best_topo && cproxy < best_cut_proxy) ||
                (topo_bad == best_topo && cproxy == best_cut_proxy &&
                 affinity > best_affinity) ||
                (topo_bad == best_topo && cproxy == best_cut_proxy &&
                 affinity == best_affinity && weight[p] < best_weight)) {
                best_topo = topo_bad;
                best_cut_proxy = cproxy;
                best_affinity = affinity;
                best_weight = weight[p];
                best_p = p;
            }
        }
        if (best_p == -1) return {};
        part[u] = best_p;
        weight[best_p] += w;
    }

    return part;
}

void build_net_side_counts_kway(const Hypergraph &hg,
                                 const vector<int> &part,
                                 int k,
                                 vector<int> &count,
                                 vector<int> &net_parts,
                                 vector<int> &net_parts_size) {
    // Flat array: count[net_idx * k + p] = pins in partition p.
    // net_parts[net_idx * k + i] = i-th occupied partition of net_idx.
    // net_parts_size[net_idx] = connectivity of net_idx (#distinct partitions).
    count.assign(hg.nets.size() * k, 0);
    net_parts.assign(hg.nets.size() * k, -1);
    net_parts_size.assign(hg.nets.size(), 0);
    for (int net_idx = 0; net_idx < static_cast<int>(hg.nets.size()); ++net_idx) {
        int base = net_idx * k;
        int &sz = net_parts_size[net_idx];
        for (int u : hg.nets[net_idx]) {
            int p = part[u];
            if (count[base + p] == 0) {
                net_parts[base + sz++] = p;
            }
            count[base + p]++;
        }
    }
}

int delta_cut_for_move_kway(const Hypergraph &hg,
                            int u, int to_part,
                            const vector<int> &part,
                            const vector<int> &count,
                            const vector<int> &net_parts_size,
                            int k) {
    int delta = 0;
    int from_part = part[u];

    for (int net_idx : hg.node_nets[u]) {
        int base = net_idx * k;
        int connectivity = net_parts_size[net_idx];
        bool before_cut = (connectivity >= 2);

        int after_conn = connectivity;
        if (count[base + from_part] == 1) after_conn--;
        if (count[base + to_part] == 0) after_conn++;
        bool after_cut = (after_conn >= 2);

        delta += (after_cut ? 1 : 0) - (before_cut ? 1 : 0);
    }
    return delta;
}

bool topo_net_bad_after_move(const TopoGraph &topo,
                             int net_idx,
                             int from_part,
                             int to_part,
                             const vector<int> &count,
                             const vector<int> &net_parts,
                             const vector<int> &net_parts_size,
                             int k) {
    int base = net_idx * k;
    int parts[128];
    int sz = 0;
    for (int i = 0; i < net_parts_size[net_idx]; ++i) {
        int p = net_parts[base + i];
        if (p == from_part && count[base + from_part] == 1) continue;
        parts[sz++] = p;
    }
    if (count[base + to_part] == 0) {
        parts[sz++] = to_part;
    }
    for (int i = 0; i < sz; ++i) {
        for (int j = i + 1; j < sz; ++j) {
            if (!topo.adj_matrix[parts[i] * topo.n + parts[j]]) {
                return true;
            }
        }
    }
    return false;
}

bool topo_net_bad_current(const TopoGraph &topo,
                          int net_idx,
                          const vector<int> &net_parts,
                          const vector<int> &net_parts_size,
                          int k) {
    int base = net_idx * k;
    int sz = net_parts_size[net_idx];
    for (int i = 0; i < sz; ++i) {
        int a = net_parts[base + i];
        for (int j = i + 1; j < sz; ++j) {
            int b = net_parts[base + j];
            if (!topo.adj_matrix[a * topo.n + b]) {
                return true;
            }
        }
    }
    return false;
}

int delta_topo_for_move_kway(const TopoGraph &topo,
                             const Hypergraph &hg,
                             int u, int to_part,
                             const vector<int> &part,
                             const vector<int> &count,
                             const vector<int> &net_parts,
                             const vector<int> &net_parts_size,
                             int k) {
    int delta = 0;
    int from_part = part[u];
    for (int net_idx : hg.node_nets[u]) {
        bool before_bad = topo_net_bad_current(topo, net_idx, net_parts,
                                               net_parts_size, k);
        bool after_bad = topo_net_bad_after_move(topo, net_idx, from_part,
                                                 to_part, count, net_parts,
                                                 net_parts_size, k);
        delta += (after_bad ? 1 : 0) - (before_bad ? 1 : 0);
    }
    return delta;
}

void apply_move_kway(const Hypergraph &hg,
                     int u, int to_part,
                     vector<int> &part,
                     vector<int> &count,
                     vector<int> &net_parts,
                     vector<int> &net_parts_size,
                     vector<int> &weight,
                     int k) {
    int from_part = part[u];
    int w = hg.node_weight[u];

    for (int net_idx : hg.node_nets[u]) {
        int base = net_idx * k;

        // Remove from_part from occupied list if it becomes empty.
        count[base + from_part]--;
        if (count[base + from_part] == 0) {
            int sz = net_parts_size[net_idx];
            for (int i = 0; i < sz; ++i) {
                if (net_parts[base + i] == from_part) {
                    net_parts[base + i] = net_parts[base + sz - 1];
                    net_parts_size[net_idx] = sz - 1;
                    break;
                }
            }
        }

        // Add to_part to occupied list if it was empty.
        if (count[base + to_part] == 0) {
            int sz = net_parts_size[net_idx];
            net_parts[base + sz] = to_part;
            net_parts_size[net_idx] = sz + 1;
        }
        count[base + to_part]++;
    }

    part[u] = to_part;
    weight[from_part] -= w;
    weight[to_part] += w;
}

// ========== k-way FM refinement ==========
//
// FM principle (k-way extension):
//   For each pass:
//     1) Build net-side count matrix count[net][p] for all k partitions
//     2) Compute initial weight vector w[p]
//     3) For each node u, compute the best gain achievable by moving u
//        to some other partition (max over k-1 targets).  Push best to heap.
//     4) Repeat up to step_limit steps:
//        - Pop best gain from heap; check balance constraints
//        - Apply move; lock the node
//        - Update gains for affected neighbors
//        - Track cumulative gain, best prefix
//     5) Roll back to best prefix; if no improvement, revert all
void fm_refine_kway(const Hypergraph &hg,
                    vector<int> &part,
                    int k,
                    mt19937 &rng,
                    int pass_limit,
                    int step_limit,
                    int candidate_size,
                    const chrono::steady_clock::time_point &deadline,
                    const TopoGraph *topo = nullptr,
                    const vector<int> *fixed_of = nullptr,
                    int topo_lambda = 0) {
    (void)rng;
    if (hg.n == 0) return;

    int total_weight = accumulate(hg.node_weight.begin(), hg.node_weight.end(), 0);
    vector<int> min_w(k), max_w(k);
    compute_balance_bounds_kway(total_weight, k, min_w, max_w);

    int steps = min(step_limit, hg.n);
    int max_heap_checks = max(64, candidate_size);
    vector<int> dedup_mark(hg.n, -1);
    int dedup_stamp = 0;
    vector<int> part_mark(k, -1);
    int part_stamp = 0;

    // 设计多pass的目的在于：第一轮 FM 所做的移动改变了图的状态，
    // 原来不可行的移动在新状态下可能可行，从而继续优化
    for (int pass = 0; pass < pass_limit; ++pass) {
        if (chrono::steady_clock::now() >= deadline) return;

        // --- Build flat net-side count array + occupied partition lists ---
        vector<int> count;
        vector<int> net_parts;
        vector<int> net_parts_size;
        build_net_side_counts_kway(hg, part, k, count, net_parts, net_parts_size);

        // --- Compute per-partition weight ---
        vector<int> weight(k, 0);
        for (int u = 0; u < hg.n; ++u) {
            weight[part[u]] += hg.node_weight[u];
        }

        vector<char> locked(hg.n, 0);
        // 固定节点锁定（初始化后再锁定，因为固定节点在 init 中已预置）
        if (fixed_of) {
            for (int u = 0; u < hg.n; ++u)
                if ((*fixed_of)[u] != -1) locked[u] = 1;
        }
        vector<MoveRecord> moves;
        moves.reserve(steps);

        // gain[u] = best gain achievable by moving u to some other partition
        // gain_target[u] = which partition gives that best gain
        vector<int> gain(hg.n, numeric_limits<int>::min());
        vector<int> gain_target(hg.n, -1);
        vector<int> gain_version(hg.n, 0);

        priority_queue<GainEntry, vector<GainEntry>, GainEntryCmp> max_heap;

        // 辅助 lambda：在 u 的相邻 partition 中找最大 gain（O(degree(u))，不依赖 k）。
        auto find_best_adjacent = [&](int u, int &out_g, int &out_t) {
            out_g = numeric_limits<int>::min();
            out_t = -1;
            part_stamp++;
            int from_part = part[u];
            auto consider = [&](int p) {
                if (p == from_part) return;
                if (p < 0 || p >= k) return;
                if (part_mark[p] == part_stamp) return;
                part_mark[p] = part_stamp;
                int delta_cut = delta_cut_for_move_kway(hg, u, p, part,
                                                        count, net_parts_size, k);
                int delta_topo = 0;
                if (topo && topo_lambda > 0) {
                    delta_topo = delta_topo_for_move_kway(*topo, hg, u, p,
                                                          part, count,
                                                          net_parts,
                                                          net_parts_size, k);
                }
                int g = -(delta_cut + topo_lambda * delta_topo);
                if (g > out_g) { out_g = g; out_t = p; }
            };
            for (int net_idx : hg.node_nets[u]) {
                int base = net_idx * k;
                for (int i = 0; i < net_parts_size[net_idx]; ++i) {
                    consider(net_parts[base + i]);
                }
            }
        };

        // 初始收益计算：只考虑与 u 有 net 相连的 partition。
        for (int u = 0; u < hg.n; ++u) {
            find_best_adjacent(u, gain[u], gain_target[u]);
            if (gain_target[u] != -1) {
                max_heap.push({gain[u], u, gain_target[u], gain_version[u]});
            }
        }

        int cumulative_gain = 0;
        int best_prefix_gain = 0;
        int best_prefix_len = 0;

        // --- Main FM loop ---
        for (int step = 0; step < steps; ++step) {
            if (chrono::steady_clock::now() >= deadline) return;

            int best_gain = numeric_limits<int>::min();
            int best_u = -1;
            int best_target = -1;

            vector<GainEntry> deferred;
            deferred.reserve(32);
            int heap_checks = 0;

            // 在最大堆中寻找最大收益候选
            while (!max_heap.empty() && heap_checks < max_heap_checks) {
                GainEntry entry = max_heap.top();
                max_heap.pop();
                heap_checks++;

                int u = entry.node;
                if (locked[u]) continue; //已锁定，跳过
                if (entry.version != gain_version[u]) continue; //版本过期，跳过

                int to_part = entry.to_part;
                int w = hg.node_weight[u];
                int from_part = part[u];

                // Check balance constraint (upper and lower bound)
                if (weight[to_part] + w > max_w[to_part]) { //目标分区超载
                    deferred.push_back(entry);
                    continue;
                }
                if (weight[from_part] - w < min_w[from_part]) { //源分区低于下限
                    deferred.push_back(entry);
                    continue;
                }

                if (topo && topo_lambda <= 0 &&
                    !topo_check_move_fast(*topo, hg, u, to_part, part,
                                                   count, net_parts, net_parts_size, k)) {
                    deferred.push_back(entry);
                    continue;
                }

                // Found a valid move
                best_u = u;
                best_gain = entry.gain;
                best_target = to_part;
                break;
            }

            for (const auto &e : deferred) max_heap.push(e); //将约束不合的部分放回堆

            if (best_u == -1) { //堆中找不到有效候选时，进行全扫描作为回退
                // Fallback: scan unlocked nodes over adjacent partitions.
                for (int u = 0; u < hg.n; ++u) {
                    if (locked[u]) continue;
                    int g, t;
                    find_best_adjacent(u, g, t);
                    if (t == -1) continue;
                    int w = hg.node_weight[u];
                    int from_part = part[u];
                    if (weight[t] + w > max_w[t]) continue;
                    if (weight[from_part] - w < min_w[from_part]) continue;
                    if (topo && topo_lambda <= 0 &&
                        !topo_check_move_fast(*topo, hg, u, t, part,
                                                       count, net_parts, net_parts_size, k)) continue;
                    if (g > best_gain) {
                        best_gain = g;
                        best_u = u;
                        best_target = t;
                    }
                }
            }

            if (best_u == -1) break;  // no valid move

            locked[best_u] = 1;
            int from_part = part[best_u];
            moves.push_back({best_u, from_part, best_target, best_gain});

            // 应用移动
            apply_move_kway(hg, best_u, best_target, part, count,
                            net_parts, net_parts_size, weight, k);

            // 更新受影响邻居增益
            dedup_stamp++;
            vector<int> affected;
            affected.reserve(64);
            //去重逻辑：best_u可能通过多个线网连接到同一个节点v
            for (int net_idx : hg.node_nets[best_u]) { //遍历移动节点连接的所有线网
                for (int v : hg.nets[net_idx]) { //遍历线网上的所有节点
                    if (locked[v]) continue;
                    if (dedup_mark[v] == dedup_stamp) continue; //去重
                    dedup_mark[v] = dedup_stamp;
                    affected.push_back(v);
                }
            }

            // 版本号机制（惰性删除）：节点v的gain更新后，堆中关于v的旧项就有了旧的版本号，需要更新
            for (int v : affected) {
                find_best_adjacent(v, gain[v], gain_target[v]);
                gain_version[v]++;
                if (gain_target[v] != -1) {
                    max_heap.push({gain[v], v, gain_target[v], gain_version[v]});
                }
            }

            // 跟踪最佳前缀，解决局部最优问题
            cumulative_gain += best_gain;
            if (cumulative_gain > best_prefix_gain) {
                best_prefix_gain = cumulative_gain;
                best_prefix_len = static_cast<int>(moves.size());
            }
        }

        if (moves.empty()) break;

        // 无改进则回滚至原点
        if (best_prefix_gain <= 0) {
            for (int i = static_cast<int>(moves.size()) - 1; i >= 0; --i) {
                auto &m = moves[i];
                // Reverse: move node back from to_part to from_part
                apply_move_kway(hg, m.node, m.from_part, part, count,
                                net_parts, net_parts_size, weight, k);
            }
            break;
        }

        // 回滚至最佳前缀
        for (int i = static_cast<int>(moves.size()) - 1; i >= best_prefix_len; --i) {
            auto &m = moves[i];
            apply_move_kway(hg, m.node, m.from_part, part, count,
                            net_parts, net_parts_size, weight, k);
        }
    }
}

struct PartitionEval {
    bool basic_valid = false;
    bool topology_valid = false;
    int cut = numeric_limits<int>::max();
    int invalid_part = 0;
    int fixed_violations = 0;
    int balance_violations = 0;
    int topo_violations = 0;
};

PartitionEval analyze_partition_kway(const Hypergraph &hg,
                                     const vector<int> &part,
                                     int k,
                                     const TopoGraph *topo,
                                     const vector<int> *fixed_of) {
    PartitionEval eval;
    eval.topology_valid = (topo == nullptr);
    if (static_cast<int>(part.size()) != hg.n) {
        eval.invalid_part = abs(static_cast<int>(part.size()) - hg.n);
        return eval;
    }

    int total_weight = accumulate(hg.node_weight.begin(), hg.node_weight.end(), 0);
    vector<int> min_w(k), max_w(k);
    compute_balance_bounds_kway(total_weight, k, min_w, max_w);
    vector<int> weight(k, 0);

    for (int u = 0; u < hg.n; ++u) {
        int p = part[u];
        if (p < 0 || p >= k) {
            eval.invalid_part++;
            continue;
        }
        weight[p] += hg.node_weight[u];
        if (fixed_of && (*fixed_of)[u] != -1 && (*fixed_of)[u] != p) {
            eval.fixed_violations++;
        }
    }

    for (int p = 0; p < k; ++p) {
        if (weight[p] < min_w[p] || weight[p] > max_w[p]) {
            eval.balance_violations++;
        }
    }

    if (eval.invalid_part || eval.fixed_violations || eval.balance_violations) {
        return eval;
    }

    eval.basic_valid = true;
    eval.cut = 0;
    vector<int> mark(k, -1);
    vector<int> parts;
    parts.reserve(k);
    int stamp = 0;
    for (const auto &net : hg.nets) {
        parts.clear();
        ++stamp;
        for (int u : net) {
            int p = part[u];
            if (mark[p] != stamp) {
                mark[p] = stamp;
                parts.push_back(p);
            }
        }
        if (parts.size() >= 2) {
            eval.cut++;
        }
        if (topo) {
            bool bad_net = false;
            for (size_t i = 0; i < parts.size() && !bad_net; ++i) {
                for (size_t j = i + 1; j < parts.size(); ++j) {
                    int a = parts[i], b = parts[j];
                    if (!topo->adj_matrix[a * topo->n + b]) {
                        eval.topo_violations++;
                        bad_net = true;
                        break;
                    }
                }
            }
        }
    }
    eval.topology_valid = (eval.topo_violations == 0);
    return eval;
}

bool validate_partition_kway(const Hypergraph &hg,
                             const vector<int> &part,
                             int k,
                             const TopoGraph *topo,
                             const vector<int> *fixed_of,
                             bool verbose = true) {
    PartitionEval eval = analyze_partition_kway(hg, part, k, topo, fixed_of);
    if (static_cast<int>(part.size()) != hg.n) {
        if (verbose) cerr << "Validation failed: partition size mismatch." << endl;
        return false;
    }
    if (!eval.basic_valid || !eval.topology_valid) {
        if (verbose) {
            cerr << "Validation failed: invalid_part=" << eval.invalid_part
                 << ", fixed_violations=" << eval.fixed_violations
                 << ", balance_violations=" << eval.balance_violations
                 << ", topo_violations=" << eval.topo_violations << endl;
        }
        return false;
    }
    return true;
}

bool partition_complete(const vector<int> &part, int n, int k) {
    if (static_cast<int>(part.size()) != n) return false;
    for (int p : part) {
        if (p < 0 || p >= k) return false;
    }
    return true;
}

struct RestartResult {
    vector<int> part;    // best partition (0..k-1) of this restart
    int cut = numeric_limits<int>::max();
    bool topology_clean = false;
};

// run_one_restart_kway -- one complete random restart for k-way partitioning.
// Coarsen -> try kInitModeCount initial partitions -> FM refine -> uncoarsen.
RestartResult run_one_restart_kway(int restart_id, //同一个id每次运行会产生完全相同的随机序列
                                    const Hypergraph &fine,
                                    int k,
                                    chrono::steady_clock::time_point deadline,
                                    const TopoGraph *topo = nullptr,
                                    const vector<int> *fixed_of = nullptr) {
    mt19937 rng(42 + restart_id * 1000);

    // --- Coarsening ---
    vector<Hypergraph> levels; //level[0]是最细粒度（原始图）
    vector<vector<int>> fine_to_coarse_maps; //[i]到[i+1]的节点映射表
    vector<vector<int>> fixed_levels;
    levels.push_back(fine);
    if (fixed_of) {
        fixed_levels.push_back(*fixed_of);
    }

    // 粗化，限制最多粗化层数以及时间阈值
    while (static_cast<int>(levels.size()) < kMaxLevels &&
           chrono::steady_clock::now() < deadline) {
        const Hypergraph &curr = levels.back();
        if (curr.n <= kMinCoarsestNodes) break;
        const vector<int> *level_fixed =
            fixed_levels.empty() ? nullptr : &fixed_levels.back();
        CoarsenResult cr = coarsen_once(curr, rng, level_fixed);
        if (cr.coarse.n <= 0 || cr.coarse.n >= curr.n) break;
        double ratio = static_cast<double>(cr.coarse.n) / static_cast<double>(curr.n); //收缩率
        fine_to_coarse_maps.push_back(std::move(cr.fine_to_coarse));
        if (fixed_of) {
            fixed_levels.push_back(std::move(cr.coarse_fixed));
        }
        levels.push_back(std::move(cr.coarse));
        if (ratio > kStopCoarsenRatio) break;
    }

    int coarsest_level = static_cast<int>(levels.size()) - 1; //最粗层索引

    const vector<int> *coarsest_fixed =
        fixed_levels.empty() ? nullptr : &fixed_levels[coarsest_level];

    // --- Try kInitModeCount initializations ---
    vector<int> local_best_part;
    int local_best_cut = numeric_limits<int>::max();
    vector<int> local_fallback_part;
    int local_fallback_cut = numeric_limits<int>::max();
    auto project_to_fine = [&](int level_idx, const vector<int> &level_part) {
        vector<int> projected = level_part;
        for (int lv = level_idx - 1; lv >= 0; --lv) {
            vector<int> finer(levels[lv].n, -1);
            const vector<int> &map = fine_to_coarse_maps[lv];
            for (int u = 0; u < levels[lv].n; ++u) {
                finer[u] = projected[map[u]];
            }
            projected.swap(finer);
        }
        return projected;
    };
    auto save_best = [&](int level_idx, const vector<int> &level_part) {
        if (!partition_complete(level_part, levels[level_idx].n, k)) return;
        vector<int> fine_part = project_to_fine(level_idx, level_part);
        if (!partition_complete(fine_part, fine.n, k)) return;
        PartitionEval eval = analyze_partition_kway(fine, fine_part, k, topo, fixed_of);
        if (!eval.basic_valid) return;
        if (!eval.topology_valid) {
            if (eval.cut < local_fallback_cut) {
                local_fallback_cut = eval.cut;
                local_fallback_part = std::move(fine_part);
            }
            return;
        }
        if (eval.cut < local_best_cut) {
            local_best_cut = eval.cut;
            local_best_part = std::move(fine_part);
        }
    };

    for (int init_mode = 0; init_mode < kInitModeCount; ++init_mode) {
        if (chrono::steady_clock::now() >= deadline) break;

        vector<int> part;
        switch (init_mode) {
        case 0:
            part = init_degree_balanced_kway(levels[coarsest_level], k, rng,
                                             topo, coarsest_fixed);
            break;
        case 1:
            part = init_random_balanced_kway(levels[coarsest_level], k, rng,
                                             topo, coarsest_fixed);
            break;
        default:
            part = init_smallnet_balanced_kway(levels[coarsest_level], k, rng,
                                               topo, coarsest_fixed);
            break;
        }

        // 在最粗层执行两轮FM优化
        // 主优化：允许大步搜索
        if (part.empty()) {
            part = init_soft_topology_balanced_kway(levels[coarsest_level],
                                                    k, rng, topo,
                                                    coarsest_fixed,
                                                    init_mode);
            if (!part.empty() &&
                !analyze_partition_kway(levels[coarsest_level], part, k,
                                        nullptr, coarsest_fixed).basic_valid) {
                part.clear();
            }
            if (part.empty()) {
                part = init_random_balanced_kway(levels[coarsest_level], k,
                                                 rng, nullptr, coarsest_fixed);
            }
        }
        if (part.empty()) continue;
        save_best(coarsest_level, part);

        fm_refine_kway(levels[coarsest_level], part, k, rng,
                       kCoarseFmPasses, kCoarseMaxSteps,
                       kCoarseCandidateSize, deadline, topo, coarsest_fixed,
                       topo ? kTopoLambdaCoarse : 0);
        save_best(coarsest_level, part);
        // 精修：小范围微调
        fm_refine_kway(levels[coarsest_level], part, k, rng,
                       kCoarsePolishPasses, kCoarsePolishMaxSteps,
                       kCoarsePolishCandidate, deadline, topo, coarsest_fixed,
                       topo ? kTopoLambdaCoarse : 0);
        save_best(coarsest_level, part);

        // 逐级细化
        for (int lv = coarsest_level - 1; lv >= 0; --lv) {
            if (chrono::steady_clock::now() >= deadline) break;
            vector<int> finer_part(levels[lv].n, 0); 
            const vector<int> &map = fine_to_coarse_maps[lv]; //细层节点在在粗层的父节点ID
            for (int u = 0; u < levels[lv].n; ++u)
                finer_part[u] = part[map[u]];//细层节点继承其父节点的分区
            part.swap(finer_part);
            const vector<int> *level_fixed =
                fixed_levels.empty() ? nullptr : &fixed_levels[lv];
            save_best(lv, part);

            auto &fm = fm_refine_kway;
            if (lv == 0) {//参数最激进，step最大
                fm(levels[lv], part, k, rng,
                   kFineFmPasses, kFineMaxSteps, kFineCandidateSize, deadline,
                   topo, level_fixed, topo ? kTopoLambdaFine : 0);
                save_best(lv, part);
                fm(levels[lv], part, k, rng,
                   kFinePolishPasses, kFinePolishMaxSteps, kFinePolishCandidate, deadline,
                   topo, level_fixed, topo ? kTopoLambdaFine : 0);
                save_best(lv, part);
            } else {
                fm(levels[lv], part, k, rng,
                   kMidFmPasses, kMidMaxSteps, kMidCandidateSize, deadline,
                   topo, level_fixed, topo ? kTopoLambdaMid : 0);
                save_best(lv, part);
                fm(levels[lv], part, k, rng,
                   kMidPolishPasses, kMidPolishMaxSteps, kMidPolishCandidate, deadline,
                   topo, level_fixed, topo ? kTopoLambdaMid : 0);
                save_best(lv, part);
            }
        }

        if (static_cast<int>(part.size()) == fine.n) save_best(0, part);
        if (chrono::steady_clock::now() >= deadline) break;
    }

    if (local_best_part.empty()) {
        vector<int> fallback =
            init_soft_topology_balanced_kway(fine, k, rng, topo, fixed_of, 1);
        if (!fallback.empty() &&
            !analyze_partition_kway(fine, fallback, k, nullptr,
                                    fixed_of).basic_valid) {
            fallback.clear();
        }
        if (fallback.empty()) {
            fallback = init_random_balanced_kway(fine, k, rng, nullptr, fixed_of);
        }
        if (!fallback.empty()) save_best(0, fallback);
    }
    if (!local_best_part.empty()) {
        return {std::move(local_best_part), local_best_cut, true};
    }
    return {std::move(local_fallback_part), local_fallback_cut, false};
}

}  // namespace

void Solution::read_benchmark(Graph &graph, string benchmark_name,
                              vector<vector<int>> *fixed_lists) {
    ifstream file(benchmark_name);

    if(!file.is_open()) {
        cerr << "Failed to open the file!" << endl;
        exit(-1);
    }

    // --- ICCAD 2021 format ---
    // Line 1: <num_nodes>
    // Line 2: <num_nets>
    // Lines 3 .. 3+num_nets-1: each line is one net, 0-based node IDs
    // Remaining lines: fixed nodes per FPGA (line i → FPGA i-1)
    int node_num, edge_num;
    string line;
    getline(file >> ws, line);
    istringstream iss(line);
    iss >> node_num;

    getline(file >> ws, line);
    istringstream iss2(line);
    iss2 >> edge_num;

    for (int node_id = 1; node_id <= node_num; ++node_id) {
        graph.get_or_create_node(node_id);
    }

    // Read nets
    for(int i = 0; i < edge_num; i++) {
        getline(file, line);
        istringstream iss(line);
        int node_id;

        Net *net = graph.add_net(i);

        while(iss >> node_id) {
            node_id++;  // 0-based → 1-based
            Node *node = graph.get_or_create_node(node_id);
            node->add_net(net);
            net->add_node(node);
        }
    }

    // Read fixed node assignments (remaining lines)
    if (fixed_lists) {
        while (getline(file, line)) {
            istringstream iss(line);
            vector<int> nodes;
            int node_id;
            while (iss >> node_id) {
                node_id++;  // 0-based → 1-based
                nodes.push_back(node_id);
            }
            fixed_lists->push_back(std::move(nodes));
        }
    }

    file.close();
}

int Solution::my_partition_algorithm(Graph graph, set<int> &X, set<int> &Y,
                                      unsigned int num_threads,
                                      int k,
                                      vector<int> *node_partition,
                                      const string &topo_file,
                                      const vector<int> *fixed_of) {
    // --- Step 1: determine parallelism ---
    if (num_threads == 0) num_threads = omp_get_max_threads();
    if (num_threads == 0) num_threads = 1;

    // --- Step 2: build read-only global data ---
    vector<int> idx_to_id;
    Hypergraph fine = build_fine_hypergraph(graph, idx_to_id);
    vector<int> fixed_internal;
    const vector<int> *fixed_internal_ptr = nullptr;
    if (fixed_of) {
        fixed_internal.assign(fine.n, -1);
        for (int u = 0; u < fine.n; ++u) {
            int node_id = idx_to_id[u];
            if (node_id >= 0 && node_id < static_cast<int>(fixed_of->size())) {
                int p = (*fixed_of)[node_id];
                if (p >= 0 && p < k) {
                    fixed_internal[u] = p;
                }
            }
        }
        fixed_internal_ptr = &fixed_internal;
    }

    // --- Step 3: adaptive time budget ---
    int64_t budget_seconds = kBudgetBaseTime;
    if (fine.n > kBudgetBaseNodes) {
        double ratio = static_cast<double>(fine.n) / static_cast<double>(kBudgetBaseNodes);
        budget_seconds = static_cast<int64_t>(static_cast<double>(kBudgetBaseTime)
                                              * pow(ratio, kBudgetExponent));
    }
    budget_seconds = max<int64_t>(budget_seconds, kBudgetMinSeconds);
    budget_seconds = min<int64_t>(budget_seconds, kBudgetMaxSeconds);
    auto deadline = chrono::steady_clock::now() + chrono::seconds(budget_seconds);

    // --- Step 3b: read topology (if provided) ---
    // TopoGraph is read-only after construction, safe for concurrent access.
    TopoGraph topo_graph;
    const TopoGraph *topo_ptr = nullptr;
    if (!topo_file.empty()) {
        topo_graph = read_topo_graph(topo_file);
        if (topo_graph.n != k) {
            cerr << "Topology node count (" << topo_graph.n
                 << ") must match k (" << k << ")." << endl;
            exit(-1);
        }
        topo_ptr = &topo_graph;
    }

    auto emit_partition = [&](const vector<int> &part) {
        X.clear();
        Y.clear();
        if (node_partition) {
            node_partition->assign(graph.get_node_num() + 1, 0);
        }
        for (int u = 0; u < fine.n; ++u) {
            int node_id = idx_to_id[u];
            int p = part[u];
            if (k == 2) {
                if (p == 0) X.insert(node_id);
                else        Y.insert(node_id);
            }
            if (node_partition) {
                if (static_cast<int>(node_partition->size()) <= node_id) {
                    node_partition->resize(node_id + 1);
                }
                (*node_partition)[node_id] = p;
            }
        }
    };

    // --- Step 4: OpenMP parallel execution ---
    // Each thread independently runs one restart via run_one_restart_kway.
    // topo_ptr and fixed_of_ptr are shared (read-only), safe.
    vector<RestartResult> results(kRestarts);

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (int r = 0; r < kRestarts; ++r) {
        results[r] = run_one_restart_kway(r, fine, k, deadline,
                                           topo_ptr, fixed_internal_ptr);
    }

    // --- Step 5: sequential merge ---
    // 顺序合并最优结果，从8个重启结果中选出cutsize最小的解
    vector<int> best_part;
    int best_cut = numeric_limits<int>::max();
    bool best_topology_clean = false;
    vector<int> fallback_part;
    int fallback_cut = numeric_limits<int>::max();
    for (auto &res : results) {
        if (res.part.empty() || res.cut == numeric_limits<int>::max()) {
            continue;
        }
        if (res.topology_clean) {
            if (res.cut >= best_cut) {
                continue;
            }
            best_cut = res.cut;
            best_part = std::move(res.part);
            best_topology_clean = true;
        } else if (!best_topology_clean && res.cut < fallback_cut) {
            fallback_cut = res.cut;
            fallback_part = std::move(res.part);
        }
    }

    // --- Step 6: fallback ---
    if (best_part.empty() && !fallback_part.empty()) {
        if (!validate_partition_kway(fine, fallback_part, k, nullptr, fixed_internal_ptr)) {
            cerr << "Fallback partition failed basic validation." << endl;
            exit(-1);
        }
        cerr << "Warning: no topology-clean partition was found before the deadline; "
             << "writing the best balance/fixed-valid fallback." << endl;
        best_cut = fallback_cut;
        best_part = std::move(fallback_part);
        best_topology_clean = false;
    }

    if (best_part.empty() || best_cut == numeric_limits<int>::max()) {
        cerr << "No partition candidate found before the deadline." << endl;
        exit(-1);
    }
    if (best_topology_clean) {
        if (!validate_partition_kway(fine, best_part, k, topo_ptr, fixed_internal_ptr)) {
            cerr << "Best partition failed final validation." << endl;
            exit(-1);
        }
    } else {
        if (!validate_partition_kway(fine, best_part, k, nullptr, fixed_internal_ptr)) {
            cerr << "Best fallback partition failed basic validation." << endl;
            exit(-1);
        }
    }

    // --- Step 7: convert result ---
    // For k=2: fill X/Y sets (backward compatible)
    // For any k: fill node_partition if requested
    emit_partition(best_part);

    return best_cut;
}
