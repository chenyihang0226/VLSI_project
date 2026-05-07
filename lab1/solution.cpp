#include "solution.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <unordered_map>
#include <vector>

namespace {
constexpr double kBalanceRatio = 0.02;
constexpr int kRestarts = 8;
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
constexpr int kTotalTimeBudgetSeconds = 90;
constexpr int kHeavyNetLimit = 64;
constexpr double kMatchDegreePenalty = 0.06;
constexpr double kSmallNetBoost = 1.7;
constexpr int kInitModeCount = 3;

struct Hypergraph {
    int n = 0;
    vector<vector<int>> nets;
    vector<vector<int>> node_nets;
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
    int min_val = static_cast<int>(ceil((0.5 - kBalanceRatio) * total_weight));
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

    vector<int> ids;
    ids.reserve(graph.get_node_num());
    for (const auto node : graph.get_nodes()) {
        ids.push_back(node->get_index());
    }
    sort(ids.begin(), ids.end());
    idx_to_id = ids;

    unordered_map<int, int> id_to_idx;
    id_to_idx.reserve(ids.size() * 2 + 1);
    for (int i = 0; i < static_cast<int>(ids.size()); ++i) {
        id_to_idx[ids[i]] = i;
    }

    hg.n = static_cast<int>(ids.size());
    hg.node_weight.assign(hg.n, 1);

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
        nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());
        if (nodes.size() >= 2) {
            hg.nets.push_back(std::move(nodes));
        }
    }

    rebuild_node_nets(hg);
    return hg;
}

CoarsenResult coarsen_once(const Hypergraph &fine, mt19937 &rng) {
    CoarsenResult result;
    int n = fine.n;
    result.fine_to_coarse.assign(n, -1);

    vector<char> unmatched(n, 1);
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);
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
            if (net_size > kHeavyNetLimit || net_size <= 1) {
                continue;
            }
            double net_weight = 1.0 / static_cast<double>(net_size - 1);
            if (net_size <= 4) {
                net_weight *= kSmallNetBoost;
            }
            for (int v : nodes) {
                if (v == u || !unmatched[v]) {
                    continue;
                }
                double deg_pen = 1.0 + kMatchDegreePenalty * static_cast<double>(fine.node_nets[v].size());
                double s = (score[v] += (net_weight / deg_pen));
                if (s > best_score) {
                    best_score = s;
                    best_v = v;
                } else if (best_v != -1 && fabs(s - best_score) < 1e-12) {
                    if (fine.node_nets[v].size() < fine.node_nets[best_v].size()) {
                        best_v = v;
                    }
                }
            }
        }

        if (best_v != -1) {
            unmatched[best_v] = 0;
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
            result.fine_to_coarse[u] = cid;
            coarse.node_weight[cid] += fine.node_weight[u];
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
        coarse_nodes.erase(unique(coarse_nodes.begin(), coarse_nodes.end()), coarse_nodes.end());
        if (coarse_nodes.size() >= 2) {
            coarse.nets.push_back(std::move(coarse_nodes));
        }
    }

    rebuild_node_nets(coarse);
    result.coarse = std::move(coarse);
    return result;
}

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
        build_net_side_counts(hg, part, count_x, count_y);

        int wx = 0;
        for (int u = 0; u < hg.n; ++u) {
            if (part[u] == 0) {
                wx += hg.node_weight[u];
            }
        }
        int wy = total_weight - wx;

        vector<char> locked(hg.n, 0);
        vector<MoveRecord> moves;
        moves.reserve(steps);

        vector<int> gain(hg.n, numeric_limits<int>::min());
        vector<int> gain_version(hg.n, 0);
        priority_queue<GainEntry, vector<GainEntry>, GainEntryCmp> max_heap;

        for (int u = 0; u < hg.n; ++u) {
            gain[u] = -delta_cut_for_move(hg, u, part, count_x, count_y);
            max_heap.push({gain[u], u, gain_version[u]});
        }

        int cumulative_gain = 0;
        int best_prefix_gain = 0;
        int best_prefix_len = 0;

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
            while (!max_heap.empty() && heap_checks < max_heap_checks) {
                GainEntry entry = max_heap.top();
                max_heap.pop();
                heap_checks++;

                int u = entry.node;
                if (locked[u]) {
                    continue;
                }
                if (entry.version != gain_version[u]) {
                    continue;
                }

                bool from_x = (part[u] == 0);
                int w = hg.node_weight[u];
                int new_wx = wx + (from_x ? -w : w);
                if (new_wx < min_x || new_wx > max_x) {
                    deferred.push_back(entry);
                    continue;
                }

                best_u = u;
                best_gain = entry.gain;
                best_from_side = part[u];
                break;
            }

            for (const auto &entry : deferred) {
                max_heap.push(entry);
            }

            if (best_u == -1) {
                // Fallback: full scan to avoid missing feasible moves due heap check cap.
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

            if (best_u == -1) {
                break;
            }

            locked[best_u] = 1;
            moves.push_back({best_u, best_from_side, best_gain});
            apply_move(hg, best_u, part, count_x, count_y, wx, wy);

            dedup_stamp++;
            vector<int> affected;
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

            for (int v : affected) {
                gain[v] = -delta_cut_for_move(hg, v, part, count_x, count_y);
                gain_version[v]++;
                max_heap.push({gain[v], v, gain_version[v]});
            }

            cumulative_gain += best_gain;
            if (cumulative_gain > best_prefix_gain) {
                best_prefix_gain = cumulative_gain;
                best_prefix_len = static_cast<int>(moves.size());
            }
        }

        if (moves.empty()) {
            break;
        }

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

}  // namespace

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

// Two-level multilevel partition:
// 1) Coarsen once
// 2) FM on coarse graph
// 3) Uncoarsen
// 4) FM on fine graph
void Solution::my_partition_algorithm(Graph graph, set<int> &X, set<int> &Y) {
    mt19937 rng(42);
    auto deadline = chrono::steady_clock::now() + chrono::seconds(kTotalTimeBudgetSeconds);

    vector<int> idx_to_id;
    Hypergraph fine = build_fine_hypergraph(graph, idx_to_id);

    vector<char> best_part;
    int best_cut = numeric_limits<int>::max();

    for (int r = 0; r < kRestarts; ++r) {
        if (chrono::steady_clock::now() >= deadline) {
            break;
        }

        vector<Hypergraph> levels;
        vector<vector<int>> fine_to_coarse_maps;
        levels.push_back(fine);

        while (static_cast<int>(levels.size()) < kMaxLevels && chrono::steady_clock::now() < deadline) {
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

        for (int init_mode = 0; init_mode < kInitModeCount; ++init_mode) {
            if (chrono::steady_clock::now() >= deadline) {
                break;
            }

            vector<char> part;
            if (init_mode == 0) {
                part = init_degree_balanced_part(levels[coarsest_level], rng);
            } else if (init_mode == 1) {
                part = init_random_balanced_part(levels[coarsest_level], rng);
            } else {
                part = init_smallnet_balanced_part(levels[coarsest_level], rng);
            }

            fm_refine_sampled(levels[coarsest_level],
                              part,
                              rng,
                              kCoarseFmPasses,
                              kCoarseMaxSteps,
                              kCoarseCandidateSize,
                              deadline);

            fm_refine_sampled(levels[coarsest_level],
                              part,
                              rng,
                              kCoarsePolishPasses,
                              kCoarsePolishMaxSteps,
                              kCoarsePolishCandidate,
                              deadline);

            for (int lv = coarsest_level - 1; lv >= 0; --lv) {
                vector<char> finer_part(levels[lv].n, 1);
                const vector<int> &map = fine_to_coarse_maps[lv];
                for (int u = 0; u < levels[lv].n; ++u) {
                    finer_part[u] = part[map[u]];
                }
                part.swap(finer_part);

                if (lv == 0) {
                    fm_refine_sampled(levels[lv],
                                      part,
                                      rng,
                                      kFineFmPasses,
                                      kFineMaxSteps,
                                      kFineCandidateSize,
                                      deadline);

                    fm_refine_sampled(levels[lv],
                                      part,
                                      rng,
                                      kFinePolishPasses,
                                      kFinePolishMaxSteps,
                                      kFinePolishCandidate,
                                      deadline);
                } else {
                    fm_refine_sampled(levels[lv],
                                      part,
                                      rng,
                                      kMidFmPasses,
                                      kMidMaxSteps,
                                      kMidCandidateSize,
                                      deadline);

                    fm_refine_sampled(levels[lv],
                                      part,
                                      rng,
                                      kMidPolishPasses,
                                      kMidPolishMaxSteps,
                                      kMidPolishCandidate,
                                      deadline);
                }

                if (chrono::steady_clock::now() >= deadline) {
                    break;
                }
            }

            int cut = evaluate_cut(fine, part);
            if (cut < best_cut) {
                best_cut = cut;
                best_part = part;
            }
        }
    }

    if (best_part.empty()) {
        best_part.assign(fine.n, 1);
    }

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
