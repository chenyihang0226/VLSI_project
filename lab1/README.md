# Lab1 Topology-Aware Hypergraph Partitioning

本目录实现了一个面向 ICCAD 2021 TopoPart benchmark 的多路超图划分程序。主程序读取 benchmark netlist、可选 FPGA topology 文件和固定节点约束，输出每个节点所属的 partition。

## 运行方式

编译：

```powershell
cd lab1
mingw32-make
```

或手动编译：

```powershell
g++ -Wall -Wextra -Werror -std=c++17 -pedantic -fopenmp -g -c solution.cpp -o solution.o
g++ -Wall -Wextra -Werror -std=c++17 -pedantic -fopenmp -g -c main.cpp -o main.o
g++ -Wall -Wextra -Werror -std=c++17 -pedantic -fopenmp -g -o main.exe main.o Graph.o Net.o Node.o solution.o evaluate.o
```

运行单个 case：

```powershell
.\main.exe "ICCAD2021-TopoPart-Benchmarks\Generated Benchmarks\case1" 8 43 --topo "ICCAD2021-TopoPart-Benchmarks\FPGA Graph\MFS2"
```

运行全部 Generated Benchmarks，使用 MFS2 topology：

```powershell
powershell -ExecutionPolicy Bypass -File .\run_mfs2_all.ps1 -Threads 8
```

脚本会生成：

- `results\caseX_partition.txt`
- `results\mfs2_logs\caseX.log`
- `results\mfs2_summary.csv`

## 当前算法结构

核心实现位于 `solution.cpp`，整体流程为：

1. 将 `Graph` 转为内部 `Hypergraph`。
2. 构建固定节点数组 `fixed_internal`。
3. 读取 topology graph，例如 `FPGA Graph\MFS2`。
4. 多线程执行若干 restart。
5. 每个 restart 内部执行：
   - 多层粗化；
   - 多种初始化；
   - k-way FM refinement；
   - 逐层 uncoarsen + refinement；
   - 保存 topology-clean 最优解，或 balance/fixed-valid fallback。
6. 全局合并所有 restart 的候选结果。
7. 输出 partition 文件。

## 已完成的实际修改

### 1. Benchmark 读取修正

`read_benchmark()` 现在会显式创建所有声明节点：

```cpp
for (int node_id = 1; node_id <= node_num; ++node_id) {
    graph.get_or_create_node(node_id);
}
```

这样即使某些节点没有出现在任何 net 中，也不会在内部图中丢失。

固定节点部分也会保留空行，保证 fixed-node line 与 FPGA partition ID 对齐。

### 2. 固定节点约束修正

原始 `fixed_of[node_id]` 会被转换为内部节点编号上的 `fixed_internal[u]`，避免原始 node id 与内部 index 混用。

粗化阶段新增 `coarse_fixed`，并禁止把 fixed 到不同 partition 的节点合并。

### 3. 平衡约束修正

k-way balance bounds 改为标准形式：

```text
[(1 - eps) * W / k, (1 + eps) * W / k]
```

避免 k 很大时某些 partition 过空。

### 4. 最终合法性验证

新增 `validate_partition_kway()` 和 `analyze_partition_kway()`，用于检查：

- partition 范围是否合法；
- fixed constraints 是否满足；
- balance 是否满足；
- topology violation 数量；
- cutsize。

其中 `analyze_partition_kway()` 将原先多次全图扫描合并为一次扫描，减少首次划分时的重复评估开销。

### 5. Anytime / fallback 结果保存

每个 restart 会在不同层级 refinement 后尝试保存候选解。

保存策略：

- topology-clean 候选进入 `local_best_part`；
- balance/fixed-valid 但 topology 不完全合法的候选进入 `local_fallback_part`；
- 全局合并时优先选择 topology-clean cut 最小解；
- 若 deadline 前没有 topology-clean 解，则输出 cut 最小的 balance/fixed-valid fallback，并打印 warning。

这避免了程序直接报：

```text
No valid partition found before the deadline.
```

### 6. FM refinement 加速结构

FM 中维护了：

- flat `count[net * k + p]`；
- `net_parts`；
- `net_parts_size`。

这样计算 cut delta 和 topology move check 时不需要每次全扫所有 partition。

FM 候选 target 主要来自节点相邻 net 已经出现过的 partition，减少了 k-way 全枚举成本。

### 7. MFS2 批量测试脚本

新增：

```text
run_mfs2_all.ps1
```

用于在当前 Windows/PowerShell 环境下批量测试 `Generated Benchmarks\case1..case8`。

脚本使用 `Start-Process` 分别重定向 stdout/stderr，避免 PowerShell 将 native stderr warning 包装成 `NativeCommandError`。

### 8. 已撤销的历史 partition 复用优化

曾短暂加入过读取已有 `results\caseX_partition.txt` 并直接复用的逻辑，用于保证 cutsize 不变并显著减少耗时。

该逻辑已经撤销。当前程序不会读取已有 partition 作为 incumbent；每次运行都会重新执行首次划分流程。

### 9. 当前存在的 soft-topology 实验性修改

当前代码中加入了两个用于进一步降低 cutsize / 改善 topology 搜索空间的实验性机制：

1. `init_soft_topology_balanced_kway()`
   - 当严格 topology 初始化失败时，使用 soft-topology 初始化；
   - 仍满足 fixed 和 capacity 上限；
   - 候选 partition 按 topology violation、局部 cut proxy、affinity、负载等因素排序；
   - 若生成结果无法通过基础合法性检查，会回退到原 random-balanced 初始化。

2. FM soft topology gain
   - 新增 `delta_topo_for_move_kway()`；
   - FM gain 从单纯 cut delta 改为：

```text
gain = -(delta_cut + lambda * delta_topo)
```

   - 当前分层 lambda 为：

```text
coarse = 1
mid    = 2
fine   = 4
```

注意：这部分是当前代码中的实际实现，但还没有完成充分质量验证。早期较大 lambda 曾导致 case1 cutsize 上升；当前较小 lambda 版本已经通过编译，但完整 case1 质量验证在运行中被中断。因此这部分后续需要继续对比：

- cutsize 是否低于原 fallback；
- topology violation 是否下降；
- 运行时间是否可接受。

## 当前注意事项

- 输出 warning `no topology-clean partition was found...` 表示最终输出的是 balance/fixed-valid fallback，不是完全 topology-clean 解。
- 当前优化重点仍在启发式搜索，结果会受到初始化、restart 和 FM 参数影响。
- 如果目标是降低首次划分耗时，不应依赖已有 `results\caseX_partition.txt`。
- 如果目标是降低 cutsize，建议优先增加诊断统计，再实现 cut-net targeted refinement 和 net-centric 主分区吸附。

