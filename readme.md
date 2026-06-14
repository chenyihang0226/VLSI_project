# VLSI Lab1 README

本文档是当前工作区的入口说明，目标是让读者快速理解：

1. 这个项目做什么；
2. 代码如何组织；
3. 如何编译和运行；
4. 当前 `solution.cpp` 实际实现了哪些算法和修改；
5. 当前结果与风险如何理解。

更长的实现细节记录见 `项目说明.md`。

## 1. 项目目标

本项目实现 ICCAD 2021 TopoPart benchmark 上的多路超图划分。

输入是一个超图：

- 节点表示电路单元；
- net 表示线网，可以连接多个节点；
- 程序将所有节点划分到 `k` 个 partition；
- 每个 partition 需要满足平衡约束；
- 优化目标是降低 `cutsize`。

cutsize 的定义：

```text
如果一条 net 的节点分布在两个或更多 partition 中，则这条 net 被切开。
cutsize = 被切开的 net 数量。
```

当指定 topology 文件时，程序还会考虑 FPGA 拓扑约束。例如 MFS2 是 43 个 FPGA 的拓扑图，因此运行 MFS2 时通常使用 `k = 43`。

## 2. 当前目录结构

```text
.
├── readme.md
├── 项目说明.md
├── 实验要求.md
├── report/
└── lab1/
    ├── main.cpp
    ├── solution.h
    ├── solution.cpp
    ├── Graph.h / Graph.cpp
    ├── Node.h  / Node.cpp
    ├── Net.h   / Net.cpp
    ├── evaluate.h / evaluate.cpp
    ├── Makefile
    ├── run_test.sh
    ├── run_mfs2_all.ps1
    ├── results/
    └── ICCAD2021-TopoPart-Benchmarks/
        ├── Generated Benchmarks/
        │   ├── case1
        │   ├── case2
        │   └── ...
        └── FPGA Graph/
            └── MFS2
```

主要文件：

| 文件 | 作用 |
| --- | --- |
| `lab1/main.cpp` | 程序入口，解析命令行，读取 benchmark，调用算法，写出结果。 |
| `lab1/solution.cpp` | 核心算法，多层次 k-way hypergraph partitioning、FM refinement、topology/fixed 处理。 |
| `lab1/solution.h` | `Solution` 接口声明。 |
| `lab1/Graph.*` | 助教框架中的图容器。 |
| `lab1/Node.*` | 节点类。 |
| `lab1/Net.*` | 超边 / 线网类。 |
| `lab1/evaluate.*` | 旧评估模块，主要用于基础 cut 评估。 |
| `lab1/run_mfs2_all.ps1` | Windows PowerShell 下批量运行 Generated Benchmarks + MFS2 的脚本。 |
| `lab1/results/` | 输出 partition 和日志。 |

## 3. 输入与输出

### 3.1 Benchmark 输入

当前主要使用 ICCAD 2021 TopoPart 格式：

```text
300000
749091
212496 250542
85805 88841
4047 221825
...
```

格式：

1. 第 1 行：节点数；
2. 第 2 行：net 数；
3. 接下来若干行：每行是一条 net，节点编号为 0-based；
4. 文件末尾：fixed-node section，每一行对应一个 FPGA / partition 的固定节点集合。

程序内部会将输入节点编号从 0-based 转为 1-based，便于和原助教框架兼容。

### 3.2 Topology 输入

例如：

```text
lab1/ICCAD2021-TopoPart-Benchmarks/FPGA Graph/MFS2
```

格式为：

```text
<num_fpga> <num_edges>
u v
u v
...
```

当前代码会检查：

```text
topology node count == k
```

因此 MFS2 需要配合 `k = 43`。

### 3.3 输出

输出文件位于：

```text
lab1/results/<case>_partition.txt
```

每一行对应一个节点：

```text
第 1 行 -> node 1 的 partition
第 2 行 -> node 2 的 partition
...
```

对于 k-way 划分，每行是 `0..k-1` 中的 partition id。

## 4. 编译与运行

### 4.1 编译

在 `lab1` 目录下：

```powershell
mingw32-make
```

也可以手动编译：

```powershell
g++ -Wall -Wextra -Werror -std=c++17 -pedantic -fopenmp -g -c main.cpp -o main.o
g++ -Wall -Wextra -Werror -std=c++17 -pedantic -fopenmp -g -c solution.cpp -o solution.o
g++ -Wall -Wextra -Werror -std=c++17 -pedantic -fopenmp -g -o main.exe main.o Graph.o Net.o Node.o solution.o evaluate.o
```

### 4.2 运行单个 case

不带 topology 的普通 k-way 划分：

```powershell
cd lab1
.\main.exe "ICCAD2021-TopoPart-Benchmarks\Generated Benchmarks\case1" 8 8
```

MFS2 topology 下的 43-way 划分：

```powershell
cd lab1
.\main.exe "ICCAD2021-TopoPart-Benchmarks\Generated Benchmarks\case1" 8 43 --topo "ICCAD2021-TopoPart-Benchmarks\FPGA Graph\MFS2"
```

参数含义：

```text
main.exe benchmark_file [num_threads] [k] [--topo topo_file]
```

- `benchmark_file`：输入 benchmark；
- `num_threads`：线程数，`0` 表示自动检测；
- `k`：partition 数，默认 2；
- `--topo`：可选 topology 文件。

### 4.3 批量运行 MFS2

在项目根目录运行：

```powershell
powershell -ExecutionPolicy Bypass -File .\lab1\run_mfs2_all.ps1 -Threads 8
```

脚本会测试：

```text
Generated Benchmarks/case1..case8
```

并固定使用：

```text
FPGA Graph/MFS2
```

输出：

```text
lab1/results/mfs2_logs/caseX.log
lab1/results/mfs2_summary.csv
lab1/results/caseX_partition.txt
```

## 5. 程序主流程

入口在 `lab1/main.cpp`。

整体流程：

```text
解析命令行参数
-> 读取 benchmark
-> 如果指定 --topo，则读取 fixed-node section
-> 构造 fixed_of[node_id]
-> 调用 Solution::my_partition_algorithm
-> 输出 results/<case>_partition.txt
-> 打印 cutsize 和耗时
```

`main.cpp` 当前不会读取已有 `results/<case>_partition.txt` 作为初始解。每次运行都会重新执行划分算法。

## 6. 核心算法概览

核心代码在 `lab1/solution.cpp`。

当前算法是多层次 k-way 超图划分：

```text
Graph
-> Hypergraph
-> coarsening
-> initial partitioning
-> FM refinement on coarse graph
-> uncoarsening
-> FM refinement on finer graphs
-> choose best result across restarts
```

### 6.1 内部 Hypergraph

内部结构：

```cpp
struct Hypergraph {
    int n;
    vector<vector<int>> nets;
    vector<vector<int>> node_nets;
    vector<int> node_weight;
};
```

含义：

- `nets[i]`：第 i 条 net 连接的内部节点编号；
- `node_nets[u]`：节点 u 连接的所有 net；
- `node_weight[u]`：节点权重，目前通常为 1。

### 6.2 Coarsening

粗化阶段将强相关节点合并为粗节点。

当前粗化会考虑 fixed-node 兼容性：

- 两个节点 fixed 到不同 partition 时禁止合并；
- fixed + unfixed 可以合并，粗节点继承 fixed partition；
- same fixed partition 可以合并。

这样避免在粗图阶段制造无法满足 fixed 约束的节点。

### 6.3 Initial Partitioning

当前有三种基础初始化：

1. degree-balanced；
2. random-balanced；
3. small-net-weighted。

在 topology 场景下，严格 topology 初始化失败时，会尝试 soft-topology 初始化。

### 6.4 FM Refinement

FM refinement 通过移动节点降低目标函数。

当前实现包括：

- balance lower/upper bound 检查；
- fixed nodes locked；
- priority queue 保存候选移动；
- lazy version 过滤堆中过期项；
- best-prefix rollback；
- affected-neighbor 局部增量更新。

为了减少 k-way 场景下的开销，FM 维护：

```text
count[net * k + p]
net_parts
net_parts_size
```

这样可以快速判断一条 net 当前跨了哪些 partition。

### 6.5 Parallel Restarts

外层 restart 使用 OpenMP 并行：

```cpp
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
for (int r = 0; r < kRestarts; ++r) {
    results[r] = run_one_restart_kway(...);
}
```

每个 restart 使用独立随机种子，结果写入 `results[r]`，没有共享写冲突。

## 7. 约束处理

### 7.1 Balance Constraint

当前使用标准 k-way balance：

```text
[(1 - eps) * W / k, (1 + eps) * W / k]
```

其中：

```text
eps = kBalanceRatio
```

这比按总权重计算 margin 更适合 k 很大的场景，例如 k=43。

### 7.2 Fixed Nodes

benchmark 末尾的 fixed-node section 会被读入。

由于输入使用原始 node id，而内部 hypergraph 使用连续编号，所以程序会构造：

```cpp
fixed_internal[u] = partition_id
```

fixed nodes 的处理：

- 初始化时预先放入指定 partition；
- FM 中 fixed nodes 被锁定；
- 粗化时 fixed 到不同 partition 的节点不能合并；
- 最终 validation 检查 fixed 约束。

### 7.3 Topology Constraint

拓扑图使用邻接矩阵表示：

```cpp
struct TopoGraph {
    int n;
    vector<bool> adj_matrix;
};
```

FM 中有快速 topology move check：

```cpp
topo_check_move_fast(...)
```

它利用 net 当前占用的 partition 集合，而不是每次扫描所有 pin 和所有 partition。

## 8. 候选结果与 Fallback

每个 restart 会保存两类候选：

```text
topology-clean candidate
balance/fixed-valid fallback candidate
```

全局合并规则：

1. 如果存在 topology-clean 解，输出 cutsize 最小的 topology-clean 解；
2. 如果 deadline 前没有 topology-clean 解，输出 cutsize 最小的 balance/fixed-valid fallback；
3. fallback 输出时会打印 warning。

典型 warning：

```text
Warning: no topology-clean partition was found before the deadline;
writing the best balance/fixed-valid fallback.
```

看到这个 warning 时，说明当前结果满足 basic constraints，但不保证完全满足 topology。

## 9. 当前已完成的实际修改

本节列出当前代码中已经实际存在的修改。

### 9.1 Benchmark 读取修正

`read_benchmark()` 显式创建所有声明节点，避免孤立节点丢失。

fixed-node section 保留空行，避免 FPGA line 和 partition id 错位。

### 9.2 fixed-node 映射修正

新增原始 node id 到内部 index 的 fixed 映射：

```text
fixed_of[node_id] -> fixed_internal[u]
```

避免 fixed 约束检查和 FM locked 逻辑用错编号。

### 9.3 粗化阶段 fixed 传播

`CoarsenResult` 新增 `coarse_fixed`。

粗化时禁止合并 fixed 到不同 partition 的节点。

### 9.4 Balance Bounds 修正

balance margin 按 `W/k` 计算，而不是按 `W` 计算。

### 9.5 一次性评估

新增：

```cpp
struct PartitionEval
PartitionEval analyze_partition_kway(...)
```

一次扫描得到：

- cutsize；
- partition id 合法性；
- fixed violations；
- balance violations；
- topology violations。

这减少了保存候选解时的重复全图扫描。

### 9.6 Anytime Candidate Saving

在粗层、各级细化后都会尝试保存当前候选。

这样即使时间预算耗尽，也有机会输出已经找到的较好结果，而不是直接无结果退出。

### 9.7 MFS2 批量测试脚本

新增：

```text
lab1/run_mfs2_all.ps1
```

该脚本会自动运行 case1 到 case8，并生成日志和汇总 CSV。

### 9.8 已撤销的历史 partition 复用

曾经加入过读取已有 `results/caseX_partition.txt` 并直接复用的逻辑。

该逻辑已经撤销，因为它只能减少已有结果情况下的重复运行时间，不能降低某个数据集第一次划分的耗时。

当前代码不会复用已有 partition。

## 10. 当前实验性修改

当前 `solution.cpp` 中存在 soft-topology 相关实验实现。

### 10.1 Soft-Topology Initialization

新增：

```cpp
init_soft_topology_balanced_kway(...)
```

当严格 topology 初始化失败时，它会尝试按以下因素选择 partition：

1. 新增 topology violation 数量；
2. 局部 cut proxy；
3. affinity；
4. partition 当前负载。

如果 soft 初始化无法通过 basic validation，则回退到 random-balanced 初始化。

### 10.2 Soft Topology FM Gain

新增：

```cpp
delta_topo_for_move_kway(...)
```

当 `topo_lambda > 0` 时，FM gain 使用：

```text
gain = -(delta_cut + topo_lambda * delta_topo)
```

当前参数：

```cpp
kTopoLambdaCoarse = 1
kTopoLambdaMid    = 2
kTopoLambdaFine   = 4
```

说明：

- 这部分代码已经存在并可编译；
- 但它仍属于实验性改动；
- 较大的 lambda 曾导致 case1 cutsize 上升；
- 当前较小 lambda 版本仍需要完整跑 case1 到 case8 验证质量。

## 11. 当前风险与注意事项

1. 如果看到 fallback warning，说明输出不是 topology-clean 解。
2. soft-topology 相关改动尚未完成充分实验验证。
3. 当前目标如果是降低首次运行耗时，应优化搜索流程本身，而不是依赖已有结果文件。
4. 当前目标如果是降低 cutsize，建议优先实现：
   - cut-net targeted refinement；
   - cut net 主分区吸附；
   - FM target partition 候选扩展；
   - topology violation 热点统计。
5. `results/` 中的 partition 文件是输出结果，不是算法默认输入。

## 12. 推荐阅读顺序

新同学建议按以下顺序阅读代码：

1. `lab1/main.cpp`：理解输入、参数和输出；
2. `lab1/solution.h`：看算法接口；
3. `lab1/solution.cpp` 中的 `read_benchmark()`；
4. `build_fine_hypergraph()`；
5. `coarsen_once()`；
6. 三个初始化函数和 `init_soft_topology_balanced_kway()`；
7. `fm_refine_kway()`；
8. `run_one_restart_kway()`；
9. `my_partition_algorithm()`。

## 13. 一句话总结

当前项目是一个多层次 k-way 超图划分器，支持 ICCAD 2021 TopoPart benchmark、fixed-node 约束、可选 FPGA topology 约束和 OpenMP 多重启并行。当前代码已经具备完整运行流程，但在 MFS2 topology-clean 结果质量和 soft-topology 策略上仍需要继续实验调优。
