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
constexpr double kBalanceRatio = 0.02;
constexpr int kRestarts = 8; //针对启发式优化，不同重启产生不同的搜索路径
constexpr int kMaxLevels = 5;
constexpr int kMinCoarsestNodes = 320;
constexpr double kStopCoarsenRatio = 0.92;
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

struct CoarsenResult {
    Hypergraph coarse;
    vector<int> fine_to_coarse;
};

struct MoveRecord {
    int node = -1;
    char from_side = 0;
    int gain = 0;
};

struct GainEntry {
    int gain = 0;
    int node = -1;
    int version = 0;
};

struct GainEntryCmp {
    bool operator()(const GainEntry &a, const GainEntry &b) const {
        return a.gain < b.gain;
    }
};

void compute_balance_bounds(int total_weight, int &min_x, int &max_x) {
    int min_val = static_cast<int>(ceil((0.5 - kBalanceRatio) * total_weight)); // Leslie：乘以权重是因为一开始原始图设定每一个节点的权重为 1， total_weight 为节点数量而粗化图的时候也是相加，所以还是相等的
    int max_val = static_cast<int>(floor((0.5 + kBalanceRatio) * total_weight));
    if (min_val < 0) {
        min_val = 0;
    }
    if (max_val > total_weight) {
        max_val = total_weight;
    }
    min_x = min_val;
    max_x = max_val;
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
CoarsenResult coarsen_once(const Hypergraph &fine, mt19937 &rng) {
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

    for (int cid = 0; cid < coarse.n; ++cid) {
        for (int u : groups[cid]) {
            result.fine_to_coarse[u] = cid; // 记录映射关系
            coarse.node_weight[cid] += fine.node_weight[u]; // 权重加和
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

// 接下来是三种划分方式
// 1.节点度排序后划分
vector<char> init_degree_balanced_part(const Hypergraph &hg, mt19937 &rng) {
    vector<char> part(hg.n, 1);
    vector<int> order(hg.n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);
    stable_sort(order.begin(), order.end(), [&hg](int a, int b) {
        return hg.node_nets[a].size() > hg.node_nets[b].size();
    });

    int total_weight = accumulate(hg.node_weight.begin(), hg.node_weight.end(), 0);
    int min_x = 0;
    int max_x = 0;
    compute_balance_bounds(total_weight, min_x, max_x);

    int wx = 0;
    int wy = 0;

    for (int u : order) {
        int w = hg.node_weight[u];
        bool put_x = false;

        if (wx + w <= min_x) {
            put_x = true;
        } else if (wy + w <= min_x) {
            put_x = false;
        } else if (wx + w > max_x) {
            put_x = false;
        } else if (wy + w > max_x) {
            put_x = true;
        } else {
            put_x = (wx <= wy);
        }

        if (put_x) {
            part[u] = 0;
            wx += w;
        } else {
            part[u] = 1;
            wy += w;
        }
    }

    return part;
}

// 2.随机划分
vector<char> init_random_balanced_part(const Hypergraph &hg, mt19937 &rng) {
    vector<char> part(hg.n, 1);
    vector<int> order(hg.n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);

    int total_weight = accumulate(hg.node_weight.begin(), hg.node_weight.end(), 0);
    int min_x = 0;
    int max_x = 0;
    compute_balance_bounds(total_weight, min_x, max_x);

    int wx = 0;
    int wy = 0;

    for (int u : order) {
        int w = hg.node_weight[u];
        bool can_x = (wx + w <= max_x);
        bool can_y = (wy + w <= max_x);
        bool put_x = false;

        if (wx < min_x && can_x) {
            put_x = true;
        } else if (wy < min_x && can_y) {
            put_x = false;
        } else if (can_x && can_y) {
            if (wx == wy) {
                put_x = ((rng() & 1U) == 0U);
            } else {
                put_x = (wx < wy);
            }
        } else if (can_x) {
            put_x = true;
        } else {
            put_x = false;
        }

        if (put_x) {
            part[u] = 0;
            wx += w;
        } else {
            part[u] = 1;
            wy += w;
        }
    }

    return part;
}

// 3.小线网保护的平衡划分
vector<char> init_smallnet_balanced_part(const Hypergraph &hg, mt19937 &rng) {
    vector<char> part(hg.n, 1);
    vector<double> score(hg.n, 0.0);

    for (const auto &net : hg.nets) {
        int net_size = static_cast<int>(net.size());
        if (net_size <= 1) {
            continue;
        }
        double w = 1.0 / static_cast<double>(net_size - 1);
        if (net_size <= 4) {
            w *= kSmallNetBoost;
        }
        for (int u : net) {
            score[u] += w;
        }
    }

    vector<int> order(hg.n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);
    stable_sort(order.begin(), order.end(), [&score](int a, int b) {
        if (fabs(score[a] - score[b]) > 1e-12) {
            return score[a] > score[b];
        }
        return a < b;
    });

    int total_weight = accumulate(hg.node_weight.begin(), hg.node_weight.end(), 0);
    int min_x = 0;
    int max_x = 0;
    compute_balance_bounds(total_weight, min_x, max_x);

    int wx = 0;
    int wy = 0;
    double sx = 0.0;
    double sy = 0.0;

    for (int u : order) {
        int w = hg.node_weight[u];
        bool can_x = (wx + w <= max_x);
        bool can_y = (wy + w <= max_x);
        bool put_x = false;

        if (wx < min_x && can_x) {
            put_x = true;
        } else if (wy < min_x && can_y) {
            put_x = false;
        } else if (can_x && can_y) {
            put_x = (sx <= sy);
        } else {
            put_x = can_x;
        }

        if (put_x) {
            part[u] = 0;
            wx += w;
            sx += score[u];
        } else {
            part[u] = 1;
            wy += w;
            sy += score[u];
        }
    }

    return part;
}

void build_net_side_counts(const Hypergraph &hg,
                           const vector<char> &part,
                           vector<int> &count_x,
                           vector<int> &count_y) {
    count_x.assign(hg.nets.size(), 0);
    count_y.assign(hg.nets.size(), 0);

    for (int net_idx = 0; net_idx < static_cast<int>(hg.nets.size()); ++net_idx) {
        for (int u : hg.nets[net_idx]) {
            if (part[u] == 0) {
                count_x[net_idx]++;
            } else {
                count_y[net_idx]++;
            }
        }
    }
}

int cut_from_counts(const vector<int> &count_x, const vector<int> &count_y) {
    int cut = 0;
    for (size_t i = 0; i < count_x.size(); ++i) {
        if (count_x[i] > 0 && count_y[i] > 0) {
            cut++;
        }
    }
    return cut;
}

int delta_cut_for_move(const Hypergraph &hg,
                       int u,
                       const vector<char> &part,
                       const vector<int> &count_x,
                       const vector<int> &count_y) {
    int delta = 0;
    bool from_x = (part[u] == 0);

    for (int net_idx : hg.node_nets[u]) {
        int cx = count_x[net_idx];
        int cy = count_y[net_idx];
        bool before_cut = (cx > 0 && cy > 0);

        if (from_x) {
            cx--;
            cy++;
        } else {
            cy--;
            cx++;
        }

        bool after_cut = (cx > 0 && cy > 0);
        delta += (after_cut ? 1 : 0) - (before_cut ? 1 : 0);
    }

    return delta;
}

void apply_move(const Hypergraph &hg,
                int u,
                vector<char> &part,
                vector<int> &count_x,
                vector<int> &count_y,
                int &wx,
                int &wy) {
    bool from_x = (part[u] == 0);
    int w = hg.node_weight[u];

    for (int net_idx : hg.node_nets[u]) {
        if (from_x) {
            count_x[net_idx]--;
            count_y[net_idx]++;
        } else {
            count_y[net_idx]--;
            count_x[net_idx]++;
        }
    }

    if (from_x) {
        part[u] = 1;
        wx -= w;
        wy += w;
    } else {
        part[u] = 0;
        wy -= w;
        wx += w;
    }
}

void rollback_move(const Hypergraph &hg,
                   int u,
                   char from_side,
                   vector<char> &part,
                   vector<int> &count_x,
                   vector<int> &count_y,
                   int &wx,
                   int &wy) {
    if (part[u] != from_side) {
        apply_move(hg, u, part, count_x, count_y, wx, wy);
    }
}

// FM原理：
// 重复若干 pass：
//     计算当前每条 net 的 X/Y 计数
//     计算每个节点移动到对侧的 gain
//     把所有节点放进最大堆

//     本 pass 内重复若干 step：
//         从堆中找 gain 最大且满足平衡约束的节点
//         移动该节点，并锁定
//         只更新受影响邻居的 gain
//         记录累计 gain 和最佳前缀

//     pass 结束：
//         如果没有正收益，全部回滚并停止
//         否则只保留最佳收益前缀，回滚后缀
void fm_refine_sampled(const Hypergraph &hg,
                       vector<char> &part,
                       mt19937 &rng,
                       int pass_limit,
                       int step_limit,
                       int candidate_size,
                       const chrono::steady_clock::time_point &deadline) {
    (void)rng;
    if (hg.n == 0) {
        return;
    }

    int total_weight = accumulate(hg.node_weight.begin(), hg.node_weight.end(), 0);
    int min_x = 0;
    int max_x = 0;
    compute_balance_bounds(total_weight, min_x, max_x);

    // 确定每轮步数与堆检查上限
    int steps = min(step_limit, hg.n);
    int max_heap_checks = max(64, candidate_size);
    vector<int> dedup_mark(hg.n, -1);
    int dedup_stamp = 0;

    for (int pass = 0; pass < pass_limit; ++pass) {
        if (chrono::steady_clock::now() >= deadline) {
            return;
        }

        vector<int> count_x;
        vector<int> count_y;
        build_net_side_counts(hg, part, count_x, count_y); // 计算出每个线网在左边和右边分别有多少节点

        // 计算出总权重 wx 和 wy
        int wx = 0;
        for (int u = 0; u < hg.n; ++u) {
            if (part[u] == 0) {
                wx += hg.node_weight[u];
            }
        }
        int wy = total_weight - wx;

        vector<char> locked(hg.n, 0);
        vector<MoveRecord> moves; //移动记录，后续方便回溯找到能收获最大gain的移动过程
        moves.reserve(steps);

        vector<int> gain(hg.n, numeric_limits<int>::min());
        vector<int> gain_version(hg.n, 0);
        
        //记录gain的最大堆
        //GainEntry是堆里的元素，成员包括gain、node索引以及版本version（为了解决maxheap不能直接更新）
        priority_queue<GainEntry, vector<GainEntry>, GainEntryCmp> max_heap; 

        for (int u = 0; u < hg.n; ++u) {
            gain[u] = -delta_cut_for_move(hg, u, part, count_x, count_y); // 如果移动到对面，能减少多少割线
            max_heap.push({gain[u], u, gain_version[u]});
        }

        int cumulative_gain = 0; //当前已经执行的移动序列的累计收益
        int best_prefix_gain = 0; //到目前为止，累计收益最大的前缀收益
        int best_prefix_len = 0; //最佳前缀包含多少步移动
        // 例如移动gain是{3，-2，4，-10}
        // 累计收益为{3.1，5，-5}
        // 则best_prefix_gain=5，best_prefix_len=3，最后保留前三步回滚第四步

        // 进行 steps 次节点移动
        for (int step = 0; step < steps; ++step) {
            if (chrono::steady_clock::now() >= deadline) {
                return;
            }

            int best_gain = numeric_limits<int>::min();
            int best_u = -1;
            char best_from_side = 0;

            vector<GainEntry> deferred;
            deferred.reserve(32);
            int heap_checks = 0;
            // 找到合适的 best_u
            while (!max_heap.empty() && heap_checks < max_heap_checks) {
                GainEntry entry = max_heap.top();
                max_heap.pop();
                heap_checks++;

                int u = entry.node;
                if (locked[u]) { // 是否已经锁定
                    continue;
                }
                if (entry.version != gain_version[u]) { // 数据版本是否过期
                    continue;
                }

                bool from_x = (part[u] == 0);
                int w = hg.node_weight[u];
                int new_wx = wx + (from_x ? -w : w);
                if (new_wx < min_x || new_wx > max_x) { // 移动后也要满足条件 min_x, max_x
                    deferred.push_back(entry); // 不满足也暂存到 数组 deferred 
                    continue;
                }

                // 由于堆顶是gain最大的，找到第一个合法候选就可以作为本步最优选择
                best_u = u;
                best_gain = entry.gain;
                best_from_side = part[u];
                break;
            }

            for (const auto &entry : deferred) {
                max_heap.push(entry);
            }

            if (best_u == -1) {
                // 执行全量扫描，以避免因堆 内存检查 上限而漏掉可行的移动步骤（前面对堆的检查是有 max_heap_checks 限制的）
                for (int u = 0; u < hg.n; ++u) {
                    if (locked[u]) {
                        continue;
                    }

                    bool from_x = (part[u] == 0);
                    int w = hg.node_weight[u];
                    int new_wx = wx + (from_x ? -w : w);
                    if (new_wx < min_x || new_wx > max_x) {
                        continue;
                    }

                    if (gain[u] > best_gain) {
                        best_gain = gain[u];
                        best_u = u;
                        best_from_side = part[u];
                    }
                }
            }

            // 没有任何合法移动就结束本轮pass
            if (best_u == -1) {
                break;
            }

            locked[best_u] = 1; // 找到合法最优后，锁定
            moves.push_back({best_u, best_from_side, best_gain}); // 记录移动清单
            apply_move(hg, best_u, part, count_x, count_y, wx, wy); // 修改它的阵营，更新体积和线网计数

            dedup_stamp++; // 时间戳
            vector<int> affected; // 受影响的邻居
            for (int net_idx : hg.node_nets[best_u]) {
                for (int v : hg.nets[net_idx]) {
                    if (locked[v]) {
                        continue;
                    }
                    if (dedup_mark[v] == dedup_stamp) {
                        continue;
                    }
                    dedup_mark[v] = dedup_stamp;
                    affected.push_back(v);
                }
            }

            // 受影响的邻居，更新 gain ，版本号 + 1，这样比重新计算所有节点的gain快很多
            for (int v : affected) {
                gain[v] = -delta_cut_for_move(hg, v, part, count_x, count_y);
                gain_version[v]++;
                max_heap.push({gain[v], v, gain_version[v]});
            }

            cumulative_gain += best_gain; // 累计收益
            if (cumulative_gain > best_prefix_gain) {
                best_prefix_gain = cumulative_gain;
                best_prefix_len = static_cast<int>(moves.size());
            }
        }

        if (moves.empty()) {
            break;
        }

        // 没有收益，回退
        if (best_prefix_gain <= 0) {
            for (int i = static_cast<int>(moves.size()) - 1; i >= 0; --i) {
                rollback_move(hg,
                              moves[i].node,
                              moves[i].from_side,
                              part,
                              count_x,
                              count_y,
                              wx,
                              wy);
            }
            break;
        }

        // 中途有一刻收益达到了巅峰，就回退到那个点
        for (int i = static_cast<int>(moves.size()) - 1; i >= best_prefix_len; --i) {
            rollback_move(hg,
                          moves[i].node,
                          moves[i].from_side,
                          part,
                          count_x,
                          count_y,
                          wx,
                          wy);
        }
    }
}

int evaluate_cut(const Hypergraph &hg, const vector<char> &part) {
    vector<int> count_x;
    vector<int> count_y;
    build_net_side_counts(hg, part, count_x, count_y);
    return cut_from_counts(count_x, count_y);
}

// ============================================================================
// Multi-thread result type
// ============================================================================

// RestartResult -- return value of a single restart.
// In the serial version, each restart's best partition was written directly
// into the outer loop's best_part/best_cut. In the multi-threaded version,
// each worker runs independently and cannot write to shared global variables.
// So we encapsulate one restart's output into this struct and return it to
// the caller (who merges results after the parallel region).
// Thread safety: each worker constructs its own RestartResult; no other
// thread accesses it concurrently, so no internal locking is needed.
struct RestartResult {
    vector<char> part;   // best partition of this restart (internal index)
    int cut = numeric_limits<int>::max();
};

// run_one_restart -- execute one complete random restart.
// Semantically equivalent to the body of the original serial loop:
//   for (int r = 0; r < kRestarts; ++r) { ... }
// Extracting it as a standalone function is the key refactoring step
// that enables multi-threading.
//
// Thread-safety design:
//   restart_id  -- by value, each thread has its own copy
//   fine        -- const ref, shared read-only; built once before any
//                  worker starts, never modified thereafter
//   deadline    -- by value (time_point is trivially copyable), each
//                  thread gets the same absolute time point
//   rng         -- thread-local mt19937, seed bound to restart_id
//   levels      -- thread-local vector, independent coarse hierarchy
//   part        -- thread-local partition vector
//
// Random seed: 42 + restart_id * 1000 guarantees:
//   1) different restart_id yields completely disjoint random sequences
//   2) deterministic reproducibility for the same restart_id
//
// Time budget: all threads share the same deadline and check it in their
// loops: if (chrono::steady_clock::now() >= deadline) break;
RestartResult run_one_restart(int restart_id,
                               const Hypergraph &fine,
                               chrono::steady_clock::time_point deadline) {
    mt19937 rng(42 + restart_id * 1000); // 每个重启使用独立随机种子，保证可复现

    // 粗化阶段
    vector<Hypergraph> levels;
    vector<vector<int>> fine_to_coarse_maps;
    levels.push_back(fine);

    while (static_cast<int>(levels.size()) < kMaxLevels &&
           chrono::steady_clock::now() < deadline) {
        const Hypergraph &curr = levels.back();
        if (curr.n <= kMinCoarsestNodes) {
            break;
        }

        CoarsenResult cr = coarsen_once(curr, rng);
        if (cr.coarse.n <= 0 || cr.coarse.n >= curr.n) {
            break;
        }

        double ratio = static_cast<double>(cr.coarse.n) / static_cast<double>(curr.n);
        fine_to_coarse_maps.push_back(std::move(cr.fine_to_coarse));
        levels.push_back(std::move(cr.coarse));

        if (ratio > kStopCoarsenRatio) {
            break;
        }
    }

    int coarsest_level = static_cast<int>(levels.size()) - 1;

    // 当前重启跟踪自己的最佳结果（无需锁）
    vector<char> local_best_part;
    int local_best_cut = numeric_limits<int>::max();

    for (int init_mode = 0; init_mode < kInitModeCount; ++init_mode) {
        if (chrono::steady_clock::now() >= deadline) {
            break;
        }

        // 初始划分
        vector<char> part;
        switch (init_mode) {
        case 0:
            part = init_degree_balanced_part(levels[coarsest_level], rng);
            break;
        case 1:
            part = init_random_balanced_part(levels[coarsest_level], rng);
            break;
        default:
            part = init_smallnet_balanced_part(levels[coarsest_level], rng);
            break;
        }

        // 粗层 FM
        fm_refine_sampled(levels[coarsest_level], part, rng,
                          kCoarseFmPasses, kCoarseMaxSteps,
                          kCoarseCandidateSize, deadline);
        fm_refine_sampled(levels[coarsest_level], part, rng,
                          kCoarsePolishPasses, kCoarsePolishMaxSteps,
                          kCoarsePolishCandidate, deadline);

        // 逐层反粗化
        for (int lv = coarsest_level - 1; lv >= 0; --lv) {
            if (chrono::steady_clock::now() >= deadline) {
                break;
            }

            vector<char> finer_part(levels[lv].n, 1);
            const vector<int> &map = fine_to_coarse_maps[lv];
            for (int u = 0; u < levels[lv].n; ++u) {
                finer_part[u] = part[map[u]];
            }
            part.swap(finer_part);

            if (lv == 0) {
                fm_refine_sampled(levels[lv], part, rng,
                                  kFineFmPasses, kFineMaxSteps,
                                  kFineCandidateSize, deadline);
                fm_refine_sampled(levels[lv], part, rng,
                                  kFinePolishPasses, kFinePolishMaxSteps,
                                  kFinePolishCandidate, deadline);
            } else {
                fm_refine_sampled(levels[lv], part, rng,
                                  kMidFmPasses, kMidMaxSteps,
                                  kMidCandidateSize, deadline);
                fm_refine_sampled(levels[lv], part, rng,
                                  kMidPolishPasses, kMidPolishMaxSteps,
                                  kMidPolishCandidate, deadline);
            }
        }

        if (chrono::steady_clock::now() >= deadline) {
            break;
        }

        int cut = evaluate_cut(fine, part);
        if (cut < local_best_cut) {
            local_best_cut = cut;
            local_best_part = std::move(part);
        }
    }

    if (local_best_part.empty()) {
        local_best_part.assign(fine.n, 1);
    }

    return {std::move(local_best_part), local_best_cut};
}

}  // namespace

void Solution::read_benchmark(Graph &graph, string benchmark_name) {
    ifstream file(benchmark_name);

    if(!file.is_open()) {
        cerr << "Failed to open the file!" << endl;
        exit(-1);
    }

    // --- ICCAD 2021 format ---
    // Line 1: <num_nodes>
    // Line 2: <num_nets>
    // Lines 3 .. 3+num_nets-1: each line is one net, <node_id> <node_id> ...
    // Remaining lines: fixed nodes per FPGA (not needed for bipartitioning)
    // Note: Node IDs are 0-based in the file; we convert to 1-based internally.
    int node_num, edge_num;
    string line;
    getline(file >> ws, line);
    istringstream iss(line);
    iss >> node_num;

    getline(file >> ws, line);
    istringstream iss2(line);
    iss2 >> edge_num;

    (void)node_num;  // informational only

    for(int i = 0; i < edge_num; i++) {
        getline(file, line);
        istringstream iss(line);
        int node_id;

        Net *net = graph.add_net(i);

        while(iss >> node_id) {
            node_id++;  // convert 0-based to 1-based indexing
            Node *node = graph.get_or_create_node(node_id);
            node->add_net(net);
            net->add_node(node);
        }
    }

    // The remaining lines (fixed nodes per FPGA) are not needed.
    // Closing the file; unread lines are harmless.
    file.close();
}

void Solution::my_partition_algorithm(Graph graph, set<int> &X, set<int> &Y,
                                      unsigned int num_threads) {
    // --- Step 1: determine parallelism ---
    // If num_threads == 0, use OpenMP's default (usually = CPU core count).
    if (num_threads == 0) {
        num_threads = omp_get_max_threads();
    }
    if (num_threads == 0) {
        num_threads = 1;
    }

    // --- Step 2: build read-only global data ---
    // fine Hypergraph: read-only after construction, shared safely among threads.
    vector<int> idx_to_id;
    Hypergraph fine = build_fine_hypergraph(graph, idx_to_id);

    // --- Step 3: adaptive time budget ---
    // Instead of a fixed wall-clock deadline, we scale the budget with the
    // problem size. Larger graphs need proportionally more time per restart.
    //
    // Formula: budget = kBudgetBaseTime * (N / kBudgetBaseNodes)^kBudgetExponent
    //   - Sub-linear exponent (0.8) accounts for coarse hierarchy reducing
    //     effective problem size at each level.
    //   - Clamped to [kBudgetMinSeconds, kBudgetMaxSeconds].
    //
    // Examples (approximate):
    //   Nodes     | Budget (s)
    //   ----------|-----------
    //    12k      |   90        (ibm01, baseline)
    //    50k      |  287        (~5 min)
    //   100k      |  499        (~8 min)
    //   300k      | 1243        (~21 min, case1 size)
    int64_t budget_seconds = kBudgetBaseTime;
    if (fine.n > kBudgetBaseNodes) {
        double ratio = static_cast<double>(fine.n) / static_cast<double>(kBudgetBaseNodes);
        budget_seconds = static_cast<int64_t>(static_cast<double>(kBudgetBaseTime)
                                              * pow(ratio, kBudgetExponent));
    }
    budget_seconds = max<int64_t>(budget_seconds, kBudgetMinSeconds);
    budget_seconds = min<int64_t>(budget_seconds, kBudgetMaxSeconds);

    auto deadline = chrono::steady_clock::now() + chrono::seconds(budget_seconds);

    // --- Step 4: OpenMP parallel execution ---
    // Pre-allocate one result slot per restart. Each thread writes to its
    // own slot (results[r]) where r is the private loop index.
    // Since different threads always access different indices, no locking
    // is needed during parallel execution.
    vector<RestartResult> results(kRestarts);

    // OpenMP parallel for:
    //   num_threads(n)  —  worker thread count
    //   schedule(dynamic)  —  finished threads pick up the next pending
    //                         restart automatically (load balancing).
    //
    // Variable sharing analysis (OpenMP default rules):
    //   results      —  shared (each thread writes results[r] with distinct r)
    //   fine         —  shared (read-only, safe)
    //   deadline     —  shared (read-only, safe)
    //   r            —  private (loop counter, handled automatically by OpenMP)
    //   run_one_restart() internals are all thread-local (stack-allocated)
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (int r = 0; r < kRestarts; ++r) {
        results[r] = run_one_restart(r, fine, deadline);
    }

    // --- Step 5: sequential merge ---
    // (Trivial compared to the parallel work; no need to parallelize.)
    vector<char> best_part;
    int best_cut = numeric_limits<int>::max();
    for (auto &res : results) {
        if (res.cut < best_cut) {
            best_cut = res.cut;
            best_part = std::move(res.part);
        }
    }

    // --- Step 6: fallback ---
    if (best_part.empty()) {
        best_part.assign(fine.n, 1);
    }

    // --- Step 7: convert result ---
    X.clear();
    Y.clear();
    for (int u = 0; u < fine.n; ++u) {
        int node_id = idx_to_id[u];
        if (best_part[u] == 0) {
            X.insert(node_id);
        } else {
            Y.insert(node_id);
        }
    }
}
